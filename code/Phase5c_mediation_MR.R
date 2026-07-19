# ============================================================
# 06c_mediation_mr.R — Phase 5 Step 2+3: 读 clump 结果, 做 MR, 算中介
# ============================================================

suppressMessages(suppressWarnings({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
  library(data.table)
}))

log_file <- "logs/06_mediation.log"
cat("\n=== Phase 5 Step 2+3: MR + 中介计算 ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n\n", file = log_file, append = TRUE)

# ---- 结局定义 ----
outcomes <- list(
  list(id = "ebi-a-GCST90018807", name = "COPD",       label = "慢阻肺"),
  list(id = "ieu-a-966",          name = "LungCancer",  label = "肺癌"),
  list(id = "ieu-a-7",            name = "CAD",         label = "冠心病")
)

# ---- 1. 读取 clump 结果 ----
cat("1. 读取 clump 结果...\n", file = log_file, append = TRUE)

step2_eqtl <- readRDS("data/step2_eqtl_all.rds")
step1_results <- readRDS("data/step1_results.rds")
sig_genes <- readRDS("data/sig_genes.rds")
cat("  基因数:", length(sig_genes), "\n", file = log_file, append = TRUE)

# 读取每个基因的 clumped 结果
clumped_data <- list()
all_clumped_snps <- character(0)

for (gene in sig_genes) {
  clumped_file <- file.path("data/clump_outputs", paste0(gene, ".clumped"))
  if (!file.exists(clumped_file)) next

  lines <- tryCatch(readLines(clumped_file), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) next

  header_idx <- grep("^\\s*CHR", lines)
  if (length(header_idx) == 0) next

  res <- tryCatch({
    read.table(text = lines[header_idx:length(lines)],
               header = TRUE, sep = "",
               stringsAsFactors = FALSE, fill = TRUE)
  }, error = function(e) NULL)

  if (is.null(res) || !"SNP" %in% names(res)) next

  clumped_snps <- res$SNP
  gene_eqtl <- step2_eqtl[GeneSymbol == gene]

  exp_dat <- data.frame(
    SNP = gene_eqtl$SNP,
    beta.exposure = as.numeric(gene_eqtl$beta_eqtl),
    se.exposure = as.numeric(gene_eqtl$se_eqtl),
    pval.exposure = as.numeric(gene_eqtl$Pvalue),
    effect_allele.exposure = gene_eqtl$AssessedAllele,
    other_allele.exposure = gene_eqtl$OtherAllele,
    eaf.exposure = as.numeric(gene_eqtl$MAF),
    samplesize.exposure = as.numeric(gene_eqtl$NrSamples),
    exposure = gene,
    units.exposure = "SD",
    id.exposure = gene,
    stringsAsFactors = FALSE
  )

  exp_clumped <- exp_dat[exp_dat$SNP %in% clumped_snps, ]
  if (nrow(exp_clumped) >= 1) {
    clumped_data[[gene]] <- exp_clumped
    all_clumped_snps <- c(all_clumped_snps, exp_clumped$SNP)
  }
}

all_clumped_snps <- unique(all_clumped_snps)
cat("  Clumped 基因数:", length(clumped_data), "\n", file = log_file, append = TRUE)
cat("  去重后 clumped SNP 数:", length(all_clumped_snps), "\n", file = log_file, append = TRUE)

if (length(clumped_data) == 0) {
  cat("\n  无 clumped 数据, 退出\n", file = log_file, append = TRUE)
  quit(save = "no")
}

# ---- 2. 提取 outcome 数据 ----
cat("\n2. 提取 outcome 数据...\n", file = log_file, append = TRUE)

safe_extract_outcome <- function(snps, id) {
  r <- NULL
  for (attempt in 1:3) {
    r <- tryCatch(extract_outcome_data(snps = snps, outcomes = id, proxies = FALSE), error = function(e) e)
    if (!inherits(r, "error") && !is.null(r) && nrow(r) > 0) return(r)
    if (inherits(r, "error")) {
      cat("    retry", attempt, ":", substr(conditionMessage(r), 1, 80), "\n", file = log_file, append = TRUE)
      Sys.sleep(5 * attempt)
    }
  }
  return(r)
}

step2_all <- list()

for (outc in outcomes) {
  cat("\n  --- 结局:", outc$label, "(", outc$id, ") ---\n", file = log_file, append = TRUE)

  out_dat <- safe_extract_outcome(all_clumped_snps, outc$id)

  if (is.null(out_dat) || inherits(out_dat, "error") || nrow(out_dat) == 0) {
    cat("    无 outcome 数据, SKIP\n", file = log_file, append = TRUE)
    next
  }
  cat("    outcome SNPs:", nrow(out_dat), "\n", file = log_file, append = TRUE)

  for (gene in names(clumped_data)) {
    exp_clumped <- clumped_data[[gene]]

    out_sub <- out_dat[out_dat$SNP %in% exp_clumped$SNP, ]
    if (nrow(out_sub) < 1) next

    dat <- tryCatch(harmonise_data(exp_clumped, out_sub, action = 2), error = function(e) {
      cat("    harmonise error for", gene, ":", substr(conditionMessage(e), 1, 60), "\n",
          file = log_file, append = TRUE)
      NULL
    })
    if (is.null(dat) || nrow(dat) == 0) next
    dat <- subset(dat, mr_keep == TRUE)
    if (nrow(dat) == 0) next

    if (nrow(dat) == 1) {
      # 单 SNP: Wald ratio = beta.outcome / beta.exposure
      b_wald <- dat$beta.outcome[1] / dat$beta.exposure[1]
      se_wald <- abs(b_wald) * sqrt((dat$se.outcome[1] / dat$beta.outcome[1])^2 +
                                     (dat$se.exposure[1] / dat$beta.exposure[1])^2)
      p_wald <- 2 * pnorm(-abs(b_wald / se_wald))

      step2_all[[length(step2_all) + 1]] <- data.frame(
        gene = gene, outcome = outc$name, outcome_label = outc$label,
        method = "Wald ratio", beta = b_wald, se = se_wald, pval = p_wald,
        nsnp = 1, note = "OK", stringsAsFactors = FALSE
      )

      cat("    ", gene, "->", outc$label, ": Wald beta=", round(b_wald, 4),
          " p=", format(p_wald, digits=2), " n=1\n",
          file = log_file, append = TRUE)
    } else {
      # 多 SNP: IVW + Egger + WME
      mr_res <- tryCatch({
        mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
      }, error = function(e) NULL)

      if (is.null(mr_res) || nrow(mr_res) == 0) next

      ivw <- mr_res[mr_res$method == "IVW", ]
      if (nrow(ivw) == 0 || length(ivw$b) == 0) next

      step2_all[[length(step2_all) + 1]] <- data.frame(
        gene = gene, outcome = outc$name, outcome_label = outc$label,
        method = "IVW", beta = ivw$b[1], se = ivw$se[1], pval = ivw$pval[1],
        nsnp = ivw$nsnp[1], note = "OK", stringsAsFactors = FALSE
      )

      for (m in c("MR Egger", "Weighted median")) {
        r <- mr_res[mr_res$method == m, ]
        if (nrow(r) > 0 && length(r$b) > 0) {
          step2_all[[length(step2_all) + 1]] <- data.frame(
            gene = gene, outcome = outc$name, outcome_label = outc$label,
            method = m, beta = r$b[1], se = r$se[1], pval = r$pval[1],
            nsnp = r$nsnp[1], note = "OK", stringsAsFactors = FALSE
          )
        }
      }

      cat("    ", gene, "->", outc$label, ": IVW beta=", round(ivw$b[1], 4),
          " p=", format(ivw$pval[1], digits=2), " n=", ivw$nsnp[1], "\n",
          file = log_file, append = TRUE)
    }
  }

  if (length(step2_all) > 0) {
    write.csv(do.call(rbind, step2_all), "results/tables/06_step2_expression_to_disease.csv", row.names = FALSE)
  }
}

if (length(step2_all) == 0) {
  cat("\n  Step 2 无结果, 退出\n", file = log_file, append = TRUE)
  quit(save = "no")
}
step2_results <- do.call(rbind, step2_all)
setDT(step2_results)
cat("\n  Step 2 完成, 总行数:", nrow(step2_results), "\n", file = log_file, append = TRUE)

# ---- 3. Step 3: 中介效应 ----
cat("\n3. Step 3: 中介效应\n", file = log_file, append = TRUE)

phase2 <- read.csv("results/tables/03_batch_mr_results.csv")
phase2_smoking <- phase2[phase2$exposure == "Smoking_initiation" & phase2$method == "IVW", c("outcome", "b")]
cat("  Phase 2 吸烟结局数:", nrow(phase2_smoking), "\n", file = log_file, append = TRUE)

mediation <- merge(
  step1_results[, .(gene = GeneSymbol, beta_step1 = beta_wald, p_step1 = p_wald)],
  step2_results[, .(gene, outcome, beta_step2 = beta, p_step2 = pval)],
  by = "gene"
)

mediation <- merge(mediation, phase2_smoking, by = "outcome", all.x = TRUE)
setnames(mediation, "b", "beta_total")

mediation[, mediation_effect := beta_step1 * beta_step2]
mediation[, prop_mediated := mediation_effect / beta_total]
mediation[, sig_mediation := p_step1 < 0.05 & p_step2 < 0.05]
mediation[, abs_prop := abs(prop_mediated)]
setorder(mediation, -abs_prop)

cat("\n  中介结果 (top 20):\n", file = log_file, append = TRUE)
for (i in 1:min(20, nrow(mediation))) {
  r <- mediation[i]
  cat("    ", r$gene, "->", r$outcome, ": med=", round(r$mediation_effect, 4),
      " prop=", round(r$prop_mediated, 3), " sig=", r$sig_mediation, "\n",
      file = log_file, append = TRUE)
}

write.csv(mediation, "results/tables/06_mediation_results.csv", row.names = FALSE)

cat("\n=== Phase 5 中介 MR 完成 ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n", file = log_file, append = TRUE)
cat("显著中介路径数:", sum(mediation$sig_mediation, na.rm = TRUE), "\n", file = log_file, append = TRUE)
