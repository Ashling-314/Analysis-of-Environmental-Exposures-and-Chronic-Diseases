# 06d_step3_only.R — 只跑 Step 3 中介效应计算
# 使用已有的 Step 2 结果 CSV
suppressMessages(suppressWarnings({
  library(data.table)
}))

log_file <- "logs/06_mediation.log"
cat("\n=== Phase 5 Step 3 (standalone): 中介效应 ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n", file = log_file, append = TRUE)

# ---- 读取 Step 1 结果 ----
step1_results <- readRDS("data/step1_results.rds")
cat("Step 1 基因数:", nrow(step1_results), "\n", file = log_file, append = TRUE)

# ---- 读取 Step 2 结果 ----
step2_results <- fread("results/tables/06_step2_expression_to_disease.csv")
cat("Step 2 行数:", nrow(step2_results), "\n", file = log_file, append = TRUE)

# ---- 读取 Phase 2 总效应 ----
phase2 <- fread("results/tables/03_batch_mr_results.csv")
phase2_smoking <- phase2[exposure == "Smoking_initiation" & method == "IVW", .(outcome, beta_total = b)]
cat("Phase 2 吸烟结局数:", nrow(phase2_smoking), "\n", file = log_file, append = TRUE)
cat("  结局列表:", paste(phase2_smoking$outcome, collapse = ", "), "\n", file = log_file, append = TRUE)

# ---- 合并 ----
mediation <- merge(
  step1_results[, .(gene = GeneSymbol, beta_step1 = beta_wald, p_step1 = p_wald)],
  step2_results[, .(gene, outcome, beta_step2 = beta, p_step2 = pval)],
  by = "gene"
)

mediation <- merge(mediation, phase2_smoking, by = "outcome", all.x = TRUE)

cat("合并后行数:", nrow(mediation), "\n", file = log_file, append = TRUE)
cat("有总效应的行数:", sum(!is.na(mediation$beta_total)), "\n", file = log_file, append = TRUE)

# ---- 计算中介效应 ----
mediation[, mediation_effect := beta_step1 * beta_step2]
mediation[, prop_mediated := mediation_effect / beta_total]
mediation[, sig_mediation := p_step1 < 0.05 & p_step2 < 0.05]
mediation[, abs_prop := abs(prop_mediated)]
setorder(mediation, -abs_prop)

# ---- 输出 top 20 ----
cat("\n  中介结果 (top 20):\n", file = log_file, append = TRUE)
for (i in 1:min(20, nrow(mediation))) {
  r <- mediation[i]
  cat("    ", r$gene, "->", r$outcome, ": med=", round(r$mediation_effect, 4),
      " prop=", round(r$prop_mediated, 3), " sig=", r$sig_mediation, "\n",
      file = log_file, append = TRUE)
}

# ---- 保存 ----
write.csv(mediation, "results/tables/06_mediation_results.csv", row.names = FALSE)

cat("\n=== Phase 5 中介 MR 完成 ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n", file = log_file, append = TRUE)
cat("总中介路径数:", nrow(mediation), "\n", file = log_file, append = TRUE)
cat("显著中介路径数:", sum(mediation$sig_mediation, na.rm = TRUE), "\n", file = log_file, append = TRUE)
cat("有总效应的中介路径数:", sum(!is.na(mediation$beta_total)), "\n", file = log_file, append = TRUE)
