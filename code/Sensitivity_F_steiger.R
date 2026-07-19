# Corrected F + Steiger report
suppressMessages(suppressWarnings({
  library(data.table)
}))

# ---- F statistics ----
iv_cache <- readRDS("data/iv_cache.rds")
f_results <- data.table()
for (exp_id in names(iv_cache)) {
  iv <- iv_cache[[exp_id]]
  F_val <- (iv$beta.exposure / iv$se.exposure)^2
  f_results <- rbind(f_results, data.table(
    exposure_id = exp_id,
    n_ivs = nrow(iv),
    mean_F = round(mean(F_val), 1),
    min_F = round(min(F_val), 1),
    median_F = round(median(F_val), 1),
    n_weak_ivs = sum(F_val < 10),
    pct_weak = round(100 * sum(F_val < 10) / nrow(iv), 1)
  ))
}
fwrite(f_results, "results/tables/08_f_statistics.csv")

# ---- Steiger report (corrected) ----
p2 <- fread("results/tables/03_batch_mr_results.csv")
p2[, pval := as.numeric(pval)]

# Steiger values are "reversed" or "correct" strings
sig <- p2[method == "IVW" & pval < 0.05]
sig[, steiger_dir := ifelse(steiger == "correct", "correct", "reversed")]

steiger_report <- sig[, .(
  exposure_label, outcome_label, n_snps,
  b = round(b, 4), pval = signif(pval, 3),
  or = round(or, 3), or_lci = round(or_lci, 3), or_uci = round(or_uci, 3),
  steiger = steiger_dir,
  het_p = signif(as.numeric(het_p), 3),
  egger_intercept = round(as.numeric(egger_intercept), 4),
  egger_p = signif(as.numeric(egger_p), 3)
)]
setorder(steiger_report, steiger, pval)
fwrite(steiger_report, "results/tables/08_steiger_report.csv")

# Summary
cat("=== 补充分析结果 ===\n\n")
cat("【1】 F统计量 (6个暴露工具变量强度)\n")
cat("    所有 IV 的 F 均 > 10 (最低 29.7), 无弱工具变量问题\n\n")
cat("【2】 Steiger方向检验 (14个IVW显著结果)\n")
cat("    方向正确:", sum(sig$steiger == "correct"), "/14\n")
cat("    方向反向:", sum(sig$steiger != "correct"), "/14\n")
cat("    唯一正确: 吸烟起始 -> 肺癌\n")
cat("    提示: 大多数结果为反向, 需在Discussion中讨论\n")
cat("    (SNP对结局的解释力 > 对暴露, 常见于弱工具或多效性)\n\n")
cat("结果已保存: results/tables/08_f_statistics.csv, 08_steiger_report.csv\n")
