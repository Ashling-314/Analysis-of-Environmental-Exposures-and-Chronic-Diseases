# ============================================================
# 05_mvmr.R — Phase 3 多变量 MR (MVMR)
# 相关暴露联合校正，分离独立因果效应
# 组合: 吸烟+饮酒, 咖啡+水果, NO2+PM2.5
# ============================================================

suppressMessages(suppressWarnings({
  library(TwoSampleMR)
  library(ieugwasr)
  library(MendelianRandomization)
  library(dplyr)
}))

log_file <- "logs/05_mvmr.log"
cat("=== Phase 3 MVMR 开始 ===\n", file = log_file)
cat("时间:", format(Sys.time()), "\n\n", file = log_file, append = TRUE)

# ---- 暴露组合定义 ----
mvmr_groups <- list(
  list(
    name = "Smoking_vs_Alcohol",
    exp1 = list(id = "ieu-b-4877",  name = "Smoking",  label = "吸烟"),
    exp2 = list(id = "ukb-b-5779",  name = "Alcohol",  label = "饮酒"),
    outcomes = c("ieu-a-7" = "CAD", "ieu-a-966" = "LungCancer",
                 "ebi-a-GCST90018807" = "COPD", "ebi-a-GCST006867" = "T2D",
                 "ebi-a-GCST009541" = "HeartFailure")
  ),
  list(
    name = "Coffee_vs_Fruit",
    exp1 = list(id = "ukb-b-5237",  name = "Coffee",  label = "咖啡"),
    exp2 = list(id = "ukb-b-3881",  name = "Fruit",   label = "水果"),
    outcomes = c("ieu-a-7" = "CAD", "ebi-a-GCST006867" = "T2D",
                 "ebi-a-GCST009541" = "HeartFailure")
  ),
  list(
    name = "NO2_vs_PM25",
    exp1 = list(id = "ukb-b-9942",   name = "NO2",   label = "NO2"),
    exp2 = list(id = "ukb-b-10817",  name = "PM25",  label = "PM2.5"),
    outcomes = c("ebi-a-GCST90018864" = "IschStroke", "ieu-a-7" = "CAD")
  )
)

# ---- 安全提取函数 ----
safe_extract_iv <- function(id) {
  r <- NULL
  for (attempt in 1:3) {
    r <- tryCatch(
      extract_instruments(outcomes = id, p1 = 5e-8, clump = TRUE),
      error = function(e) e
    )
    if (!inherits(r, "error") && !is.null(r) && nrow(r) > 0) {
      return(r)
    }
    if (inherits(r, "error")) {
      msg <- conditionMessage(r)
      cat("  retry", attempt, ":", substr(msg, 1, 50), "\n", file = log_file, append = TRUE)
      if (grepl("401", msg)) Sys.sleep(5 * attempt)
    }
    Sys.sleep(2)
  }
  return(r)
}

safe_extract_outcome <- function(snps, id) {
  r <- NULL
  for (attempt in 1:3) {
    r <- tryCatch(
      extract_outcome_data(snps = snps, outcomes = id),
      error = function(e) e
    )
    if (!inherits(r, "error") && !is.null(r) && nrow(r) > 0) {
      return(r)
    }
    if (inherits(r, "error")) {
      msg <- conditionMessage(r)
      cat("  retry", attempt, ":", substr(msg, 1, 50), "\n", file = log_file, append = TRUE)
      if (grepl("401", msg)) Sys.sleep(5 * attempt)
    }
    Sys.sleep(2)
  }
  return(r)
}

# ---- MVMR 主函数 ----
run_mvmr <- function(exp1_id, exp2_id, outcome_id, outcome_name) {
  cat("\n--- MVMR:", exp1_id, "+", exp2_id, "->", outcome_id, "(", outcome_name, ") ---\n",
      file = log_file, append = TRUE)

  # 错误返回的统一格式
  err_df <- function(note, n = 0) data.frame(
    outcome = outcome_name, method = "MVMR-IVW",
    exp1_beta = NA, exp1_se = NA, exp1_p = NA,
    exp2_beta = NA, exp2_se = NA, exp2_p = NA,
    n_snps = n, note = note)

  # 1. 提取各暴露 IV
  exp1_iv <- safe_extract_iv(exp1_id)
  exp2_iv <- safe_extract_iv(exp2_id)

  if (is.null(exp1_iv) || nrow(exp1_iv) == 0) {
    cat("  exp1 IVs: 0, SKIP\n", file = log_file, append = TRUE)
    return(err_df("no_exp1_IVs"))
  }
  if (is.null(exp2_iv) || nrow(exp2_iv) == 0) {
    cat("  exp2 IVs: 0, SKIP\n", file = log_file, append = TRUE)
    return(err_df("no_exp2_IVs"))
  }

  cat("  exp1 IVs:", nrow(exp1_iv), " exp2 IVs:", nrow(exp2_iv), "\n",
      file = log_file, append = TRUE)

  # 2. 合并 SNP 列表（去重）
  all_snps <- unique(c(exp1_iv$SNP, exp2_iv$SNP))
  cat("  combined SNPs:", length(all_snps), "\n", file = log_file, append = TRUE)

  # 3. 提取各暴露对所有 SNP 的效应量（extract_outcome_data 返回 .outcome 列名，需重命名）
  rename_to_exposure <- function(dat, exp_name, exp_id) {
    if (is.null(dat) || nrow(dat) == 0) return(NULL)
    # 重命名 .outcome -> .exposure
    names(dat) <- gsub("\\.outcome$", ".exposure", names(dat))
    dat$id.exposure <- exp_id
    dat$exposure <- exp_name
    dat$units.exposure <- "SD"
    return(dat)
  }

  exp1_full <- rename_to_exposure(safe_extract_outcome(all_snps, exp1_id), "exp1", exp1_id)
  exp2_full <- rename_to_exposure(safe_extract_outcome(all_snps, exp2_id), "exp2", exp2_id)
  out_dat   <- safe_extract_outcome(all_snps, outcome_id)

  if (is.null(exp1_full) || is.null(exp2_full) || is.null(out_dat)) {
    cat("  data extraction failed, SKIP\n", file = log_file, append = TRUE)
    return(err_df("extraction_failed"))
  }

  cat("  exp1_full:", nrow(exp1_full), " exp2_full:", nrow(exp2_full),
      " out_dat:", nrow(out_dat), "\n", file = log_file, append = TRUE)

  # 4. Harmonise 各暴露与结局
  dat1 <- harmonise_data(exp1_full, out_dat, action = 2)
  dat2 <- harmonise_data(exp2_full, out_dat, action = 2)

  dat1 <- subset(dat1, mr_keep == TRUE)
  dat2 <- subset(dat2, mr_keep == TRUE)

  # 5. 取三方面都有的 SNP（交集）
  common_snps <- intersect(dat1$SNP, dat2$SNP)
  if (length(common_snps) < 3) {
    cat("  common SNPs:", length(common_snps), " < 3, SKIP\n", file = log_file, append = TRUE)
    return(err_df("insufficient_common_SNPs", length(common_snps)))
  }

  dat1 <- dat1[dat1$SNP %in% common_snps, ]
  dat2 <- dat2[dat2$SNP %in% common_snps, ]

  # 按 SNP 排序保证对齐
  dat1 <- dat1[order(dat1$SNP), ]
  dat2 <- dat2[order(dat2$SNP), ]

  cat("  common SNPs after harmonise:", length(common_snps), "\n",
      file = log_file, append = TRUE)

  # 6. 构造 MVMR 输入并运行
  mvmv_input <- tryCatch({
    mr_mvinput(
      bx = cbind(dat1$beta.exposure, dat2$beta.exposure),
      bxse = cbind(dat1$se.exposure, dat2$se.exposure),
      by = dat1$beta.outcome,
      byse = dat1$se.outcome
    )
  }, error = function(e) {
    cat("  mvinput error:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    return(NULL)
  })

  if (is.null(mvmv_input)) {
    return(err_df("mvinput_failed", length(common_snps)))
  }

  # IVW
  ivw_res <- tryCatch(mr_mvivw(mvmv_input), error = function(e) {
    cat("  mvivw error:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    return(NULL)
  })

  # Egger (如果 SNP 足够)
  egger_res <- NULL
  if (length(common_snps) >= 5) {
    egger_res <- tryCatch(mr_mvegger(mvmv_input), error = function(e) NULL)
  }

  # 7. 整理结果
  results <- list()

  if (!is.null(ivw_res)) {
    results[["IVW"]] <- data.frame(
      outcome = outcome_name,
      method = "MVMR-IVW",
      exp1_beta = ivw_res@Estimate[1],
      exp1_se = ivw_res@StdError[1],
      exp1_p = ivw_res@Pvalue[1],
      exp2_beta = ivw_res@Estimate[2],
      exp2_se = ivw_res@StdError[2],
      exp2_p = ivw_res@Pvalue[2],
      n_snps = length(common_snps),
      note = "OK"
    )
    cat("  IVW: exp1 b=", round(ivw_res@Estimate[1], 4),
        " p=", format(ivw_res@Pvalue[1], digits = 3),
        " | exp2 b=", round(ivw_res@Estimate[2], 4),
        " p=", format(ivw_res@Pvalue[2], digits = 3), "\n",
        file = log_file, append = TRUE)
  }

  if (!is.null(egger_res)) {
    # MVEgger 的 slot 名可能不同，用 tryCatch 安全提取
    egger_df <- tryCatch({
      sn <- slotNames(egger_res)
      est_slot <- if ("Estimate" %in% sn) egger_res@Estimate else rep(NA, 3)
      se_slot <- if ("StdError" %in% sn) egger_res@StdError else
                 if ("StdErr" %in% sn) egger_res@StdErr else rep(NA, 3)
      pv_slot <- if ("Pvalue" %in% sn) egger_res@Pvalue else
                 if ("Pval" %in% sn) egger_res@Pval else rep(NA, 3)
      # Egger: 第一个是截距, 后两个是暴露效应
      data.frame(
        outcome = outcome_name,
        method = "MVMR-Egger",
        exp1_beta = est_slot[2],
        exp1_se = se_slot[2],
        exp1_p = pv_slot[2],
        exp2_beta = est_slot[3],
        exp2_se = se_slot[3],
        exp2_p = pv_slot[3],
        n_snps = length(common_snps),
        note = "OK"
      )
    }, error = function(e) {
      cat("  Egger extract error:", conditionMessage(e), "\n", file = log_file, append = TRUE)
      return(NULL)
    })
    if (!is.null(egger_df)) results[["Egger"]] <- egger_df
  }

  if (length(results) == 0) {
    return(data.frame(outcome = outcome_name, method = "MVMR-IVW",
                      exp1_beta = NA, exp1_se = NA, exp1_p = NA,
                      exp2_beta = NA, exp2_se = NA, exp2_p = NA,
                      n_snps = length(common_snps), note = "all_methods_failed"))
  }

  do.call(rbind, results)
}

# ---- 主循环 ----
all_mvmr <- list()
cat("\n开始 MVMR 分析...\n", file = log_file, append = TRUE)

for (grp in mvmr_groups) {
  cat("\n========================================\n", file = log_file, append = TRUE)
  cat("Group:", grp$name, "\n", file = log_file, append = TRUE)
  cat("========================================\n", file = log_file, append = TRUE)

  for (outc_id in names(grp$outcomes)) {
    outc_name <- grp$outcomes[outc_id]
    res <- run_mvmr(grp$exp1$id, grp$exp2$id, outc_id, outc_name)
    res$group <- grp$name
    res$exp1_name <- grp$exp1$label
    res$exp2_name <- grp$exp2$label
    all_mvmr[[length(all_mvmr) + 1]] <- res

    # 增量保存
    combined <- do.call(rbind, all_mvmr)
    write.csv(combined, "results/tables/05_mvmr_results.csv", row.names = FALSE)
    cat("  -> 已保存到 CSV\n", file = log_file, append = TRUE)
  }
}

cat("\n=== Phase 3 MVMR 完成 ===\n", file = log_file, append = TRUE)
cat("时间:", format(Sys.time()), "\n", file = log_file, append = TRUE)
cat("总结果行数:", nrow(do.call(rbind, all_mvmr)), "\n", file = log_file, append = TRUE)
