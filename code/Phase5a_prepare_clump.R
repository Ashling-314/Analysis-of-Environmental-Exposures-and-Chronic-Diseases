# ============================================================
# 06a_prepare_clump.R — Phase 5 Step 1: 准备 LD clump 输入
# 读取 eQTL 数据, 计算步骤1 Wald ratio, 写出每个基因的 SNP 列表
# ============================================================

suppressMessages(suppressWarnings({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
  library(data.table)
}))

log_file <- "logs/06_mediation.log"
cat("=== Phase 5 Step 1: 准备 clump 输入 ===\n", file = log_file)
cat("时间:", format(Sys.time()), "\n\n", file = log_file, append = TRUE)

# ---- 1. 读取数据 ----
cat("1. 读取数据...\n", file = log_file, append = TRUE)

smoking_iv <- readRDS("data/smoking_iv_full.rds")
setDT(smoking_iv)
cat("  吸烟 IVs:", nrow(smoking_iv), "\n", file = log_file, append = TRUE)

eqtl_hits <- fread("data/eqtl_smoking_hits.tsv", header = FALSE)
setnames(eqtl_hits, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele",
                       "OtherAllele","Zscore","Gene","GeneSymbol","GeneChr",
                       "GenePos","NrCohorts","NrSamples","FDR","BonferroniP"))
cat("  吸烟-eQTL 重叠:", nrow(eqtl_hits), " 条\n", file = log_file, append = TRUE)

eqtl_full <- fread("data/eqtl_mediator_full.tsv", header = FALSE)
setnames(eqtl_full, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele",
                       "OtherAllele","Zscore","Gene","GeneSymbol","GeneChr",
                       "GenePos","NrCohorts","NrSamples","FDR","BonferroniP"))
cat("  全部 cis-eQTL:", nrow(eqtl_full), " 行\n", file = log_file, append = TRUE)

af_data <- fread("data/mediator_af.tsv", header = FALSE)
setnames(af_data, c("SNP","chr","pos","AlleleA","AlleleB",
                     "allA","allAB","allB","AlleleB_freq"))
af_data[, MAF := pmin(AlleleB_freq, 1 - AlleleB_freq)]

# ---- 2. 计算 eQTL beta/se ----
cat("\n2. 计算 eQTL beta/se...\n", file = log_file, append = TRUE)

eqtl_hits <- merge(eqtl_hits, af_data[, .(SNP, MAF)], by = "SNP", all.x = TRUE)
eqtl_full <- merge(eqtl_full, af_data[, .(SNP, MAF)], by = "SNP", all.x = TRUE)
eqtl_hits[is.na(MAF), MAF := 0.3]
eqtl_full[is.na(MAF), MAF := 0.3]

eqtl_hits[, beta_eqtl := Zscore / sqrt(NrSamples * 2 * MAF * (1 - MAF))]
eqtl_hits[, se_eqtl := 1 / sqrt(NrSamples * 2 * MAF * (1 - MAF))]
eqtl_full[, beta_eqtl := Zscore / sqrt(NrSamples * 2 * MAF * (1 - MAF))]
eqtl_full[, se_eqtl := 1 / sqrt(NrSamples * 2 * MAF * (1 - MAF))]

# ---- 3. Step 1: 吸烟 -> 基因表达 (Wald ratio) ----
cat("\n3. Step 1: 吸烟 -> 基因表达\n", file = log_file, append = TRUE)

step1_data <- merge(
  eqtl_hits[, .(SNP, GeneSymbol, beta_eqtl, se_eqtl, Pvalue, AssessedAllele, OtherAllele)],
  smoking_iv[, .(SNP, beta.exposure, se.exposure, pval.exposure, effect_allele.exposure, other_allele.exposure)],
  by = "SNP"
)

step1_data[, flip := !(toupper(AssessedAllele) == toupper(effect_allele.exposure) &
                        toupper(OtherAllele) == toupper(other_allele.exposure))]
step1_data[, flip2 := !(toupper(AssessedAllele) == toupper(other_allele.exposure) &
                         toupper(OtherAllele) == toupper(effect_allele.exposure))]
step1_data[, aligned := flip | flip2]
step1_data[, need_flip := flip & !flip2]
step1_data[need_flip == TRUE, beta_eqtl := -beta_eqtl]
step1_data <- step1_data[aligned == TRUE]

step1_data[, beta_wald := beta_eqtl / beta.exposure]
step1_data[, se_wald := abs(beta_wald) * sqrt((se_eqtl/beta_eqtl)^2 + (se.exposure/beta.exposure)^2)]
step1_data[, z_wald := beta_wald / se_wald]
step1_data[, p_wald := 2 * pnorm(-abs(z_wald))]

step1_results <- step1_data[, .(GeneSymbol, SNP, beta_wald, se_wald, p_wald,
                                 beta_smoking = beta.exposure, beta_eqtl)]
setorder(step1_results, p_wald)

sig_genes <- step1_results[p_wald < 0.05, unique(GeneSymbol)]
cat("  Step 1 显著基因:", length(sig_genes), "\n", file = log_file, append = TRUE)

write.csv(step1_results, "results/tables/06_step1_smoking_to_expression.csv", row.names = FALSE)

if (length(sig_genes) < 5) {
  sig_genes <- head(step1_results[, unique(GeneSymbol)], 20)
}

# ---- 4. 写出每个基因的 SNP 列表 + exposure 数据 ----
cat("\n4. 写出 clump 输入文件...\n", file = log_file, append = TRUE)

step2_eqtl <- eqtl_full[GeneSymbol %in% sig_genes]

# 创建 clump 输入目录
clump_dir <- "data/clump_inputs"
dir.create(clump_dir, showWarnings = FALSE)

# 清理旧文件
old_files <- list.files(clump_dir, full.names = TRUE)
if (length(old_files) > 0) file.remove(old_files)

gene_list <- data.frame(gene = character(), n_snps = integer(), stringsAsFactors = FALSE)

for (gene in sig_genes) {
  gene_eqtl <- step2_eqtl[GeneSymbol == gene]

  # 写 PLINK clump 输入文件
  clump_file <- file.path(clump_dir, paste0(gene, ".txt"))
  write.table(data.frame(SNP = gene_eqtl$SNP, P = gene_eqtl$Pvalue),
              file = clump_file, sep = "\t", quote = FALSE, row.names = FALSE)

  gene_list <- rbind(gene_list, data.frame(gene = gene, n_snps = nrow(gene_eqtl), stringsAsFactors = FALSE))
}

# 保存基因列表和 exposure 数据
write.csv(gene_list, "data/clump_gene_list.csv", row.names = FALSE)

# 保存所有基因的 exposure 数据（clump 后用于 MR）
saveRDS(step2_eqtl, "data/step2_eqtl_all.rds")
saveRDS(step1_results, "data/step1_results.rds")
saveRDS(sig_genes, "data/sig_genes.rds")

cat("  写出", length(sig_genes), "个基因的 clump 输入文件到", clump_dir, "\n", file = log_file, append = TRUE)
cat("  总 SNP 数:", sum(gene_list$n_snps), "\n", file = log_file, append = TRUE)

cat("\n=== Step 1 完成, 请运行 06b_run_clump.sh ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n", file = log_file, append = TRUE)
