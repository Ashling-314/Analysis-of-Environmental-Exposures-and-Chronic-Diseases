# Comprehensive sensitivity analysis report from existing batch_mr_results
suppressMessages(suppressWarnings({
  library(data.table)
}))

p2 <- fread("results/tables/03_batch_mr_results.csv")
p2[, pval := as.numeric(pval)]
p2[, presso_global_p := as.numeric(presso_global_p)]
p2[, het_p := as.numeric(het_p)]
p2[, egger_intercept := as.numeric(egger_intercept)]
p2[, egger_p := as.numeric(egger_p)]

sig <- p2[method == "IVW" & pval < 0.05]
sig[, steiger_dir := ifelse(steiger == "correct", "correct", "reversed")]

# Comprehensive sensitivity table
sens_report <- sig[, .(
  exposure_label, outcome_label, n_snps,
  ivw_beta = round(b, 4),
  ivw_p = signif(pval, 3),
  ivw_or = round(or, 3),
  min_F = fstat,
  # Steiger
  steiger = steiger_dir,
  # Heterogeneity
  Q = round(het_q, 2),
  Q_p = signif(het_p, 3),
  het_significant = het_p < 0.05,
  # MR-Egger
  egger_intercept = round(egger_intercept, 4),
  egger_p = signif(egger_p, 3),
  egger_suggests_pleiotropy = egger_p < 0.05,
  # MR-PRESSO
  presso_global_p = signif(presso_global_p, 3),
  presso_outlier = presso_outlier,
  presso_significant = presso_global_p < 0.05
)]
setorder(sens_report, ivw_p)
fwrite(sens_report, "results/tables/08_sensitivity_report.csv")

cat("=== 敏感性分析综合报告 (14个IVW显著结果) ===\n\n")
cat("异质性 (Q_p < 0.05):", sum(sens_report$het_significant), "/14\n")
cat("MR-Egger 截距显著 (提示水平多效性):", sum(sens_report$egger_suggests_pleiotropy), "/14\n")
cat("MR-PRESSO 全局检验显著:", sum(sens_report$presso_significant), "/14\n")
cat("MR-PRESSO 检测到离群点:", sum(sens_report$presso_outlier, na.rm = TRUE), "/14\n")
cat("Steiger 方向正确:", sum(sens_report$steiger == "correct"), "/14\n\n")

cat("--- 异质性/Egger/PRESSO 均显著的对 (需特别注意) ---\n")
problematic <- sens_report[het_significant == TRUE & presso_significant == TRUE]
if (nrow(problematic) > 0) {
  for (i in 1:nrow(problematic)) {
    r <- problematic[i]
    cat(sprintf("  %s -> %s: Q_p=%.3f, PRESSO_p=%.3f, outlier=%s\n",
                r$exposure_label, r$outcome_label, r$Q_p, r$presso_global_p, r$presso_outlier))
  }
} else {
  cat("  无\n")
}

cat("\n结果已保存: results/tables/08_sensitivity_report.csv\n")
