# ============================================================
# 03_batch_mr.R v2 — 批量单变量 MR（Phase 2）
# 6 暴露 x 8 结局 = 48 对分析
# v2 修复: 去掉 return()、逐对保存、IV 缓存、401 重试
# ============================================================

suppressMessages(suppressWarnings({
  library(TwoSampleMR)
  library(ieugwasr)
  library(MendelianRandomization)
  library(MRPRESSO)
  library(dplyr)
}))

log_file <- "logs/03_batch_mr.log"
cat("=== 批量单变量 MR v2 开始 ===\n", file = log_file)
cat("时间:", format(Sys.time()), "\n\n", file = log_file, append = TRUE)

# ---- 暴露定义 ----
exposures <- list(
  list(id = "ieu-b-4877",    name = "Smoking_initiation",    label = "吸烟起始",     type = "binary"),
  list(id = "ukb-b-5779",    name = "Alcohol_intake_freq",   label = "饮酒频率",     type = "continuous"),
  list(id = "ukb-b-9942",    name = "NO2_air_pollution",     label = "NO2空气污染",  type = "continuous"),
  list(id = "ukb-b-10817",   name = "PM25_air_pollution",    label = "PM2.5空气污染",type = "continuous"),
  list(id = "ukb-b-5237",    name = "Coffee_intake",         label = "咖啡摄入",     type = "continuous"),
  list(id = "ukb-b-3881",    name = "Fresh_fruit_intake",    label = "新鲜水果摄入", type = "continuous")
)

# ---- 结局定义 ----
outcomes <- list(
  list(id = "ieu-a-7",              name = "CAD",           label = "冠心病"),
  list(id = "ieu-a-966",            name = "LungCancer",    label = "肺癌"),
  list(id = "ebi-a-GCST90018864",   name = "IschStroke",    label = "缺血性脑卒中"),
  list(id = "ebi-a-GCST009541",     name = "HeartFailure",  label = "心力衰竭"),
  list(id = "ebi-a-GCST90018807",   name = "COPD",          label = "慢阻肺"),
  list(id = "ebi-a-GCST006867",     name = "T2D",           label = "2型糖尿病"),
  list(id = "ieu-a-1187",           name = "MDD",           label = "抑郁症"),
  list(id = "ieu-b-5067",           name = "Alzheimer",     label = "阿尔茨海默病")
)

# ---- F 统计量 ----
calc_fstat <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) NA
  else mean((dat$beta.exposure / dat$se.exposure)^2, na.rm = TRUE)
}

# ---- 带 401 重试的 API 调用 ----
safe_extract_instruments <- function(outcome_id, p1 = 5e-8) {
  r <- NULL
  for (attempt in 1:3) {
    r <- tryCatch(
      extract_instruments(outcomes = outcome_id, p1 = p1, clump = TRUE),
      error = function(e) e
    )
    if (!inherits(r, "error")) return(r)
    else {
      msg <- conditionMessage(r)
      cat("    retry", attempt, ":", substr(msg, 1, 40), "\n", file = log_file, append = TRUE)
      if (grepl("401", msg)) Sys.sleep(5 * attempt)
      else break
    }
  }
  r  # 返回最后一次结果（可能是 error 对象）
}

safe_extract_outcome <- function(snps, outcome_id) {
  r <- NULL
  for (attempt in 1:3) {
    r <- tryCatch(
      extract_outcome_data(snps = snps, outcomes = outcome_id),
      error = function(e) e
    )
    if (!inherits(r, "error")) return(r)
    else {
      msg <- conditionMessage(r)
      cat("    retry", attempt, ":", substr(msg, 1, 40), "\n", file = log_file, append = TRUE)
      if (grepl("401", msg)) Sys.sleep(5 * attempt)
      else break
    }
  }
  r
}

# ---- 单对 MR 分析（纯函数，不用 return） ----
run_mr_pair <- function(exp, outc, exp_dat) {
  # exp_dat 已提前缓存，这里只提取结局
  out_dat <- safe_extract_outcome(exp_dat$SNP, outc$id)

  if (inherits(out_dat, "error") || is.null(out_dat) || nrow(out_dat) == 0) {
    df_error(exp, outc, n_iv = nrow(exp_dat), note = "no_outcome_data")
  } else {
    dat <- harmonise_data(exp_dat, out_dat, action = 2)
    dat <- subset(dat, mr_keep == TRUE)

    if (nrow(dat) < 2) {
      df_error(exp, outc, n_iv = nrow(exp_dat), note = "insufficient_after_harmonise")
    } else {
      analyze_pair(exp, outc, dat)
    }
  }
}

# ---- 实际分析逻辑 ----
analyze_pair <- function(exp, outc, dat) {
  n_iv <- nrow(dat)
  fstat <- calc_fstat(dat)

  # MR
  mr_res <- tryCatch(
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                             "mr_weighted_median", "mr_weighted_mode")),
    error = function(e) NULL
  )

  if (is.null(mr_res) || nrow(mr_res) == 0) {
    df_error(exp, outc, n_iv = n_iv, note = "mr_failed")
  } else {
    # Steiger
    steiger <- tryCatch({
      s <- directionality_test(dat)
      ifelse(s$steiger_pval < 0.05, "reversed", "correct")
    }, error = function(e) "NA")

    # 异质性 & 多效性
    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)

    # MR-PRESSO
    presso_global <- NA
    presso_outlier <- FALSE
    if (n_iv >= 4) {
      presso <- tryCatch(
        mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                  SdOutcome = "se.outcome", SdExposure = "se.exposure",
                  data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                  NbDistribution = 1000, SignifThreshold = 0.05),
        error = function(e) NULL
      )
      if (!is.null(presso)) {
        presso_global <- presso$`MR-PRESSO results`$`Global Test`$Pvalue[1]
        presso_outlier <- !is.null(presso$`MR-PRESSO results`$`Outlier Test`)
      }
    }

    # 整理 IVW 行
    ivw_row <- mr_res[mr_res$method == "Inverse variance weighted", ]
    if (nrow(ivw_row) == 0) ivw_row <- mr_res[1, ]

    b <- ivw_row$b
    se <- ivw_row$se
    pval <- ivw_row$pval

    # 主结果行
    res_df <- data.frame(
      exposure = exp$name, outcome = outc$name,
      exposure_label = exp$label, outcome_label = outc$label,
      n_snps = n_iv, method = "IVW",
      b = b, se = se, pval = pval,
      or = exp(b), or_lci = exp(b - 1.96 * se), or_uci = exp(b + 1.96 * se),
      fstat = round(fstat, 1), steiger = steiger,
      het_q = ifelse(!is.null(het) && nrow(het) > 0, het$Q[1], NA),
      het_p = ifelse(!is.null(het) && nrow(het) > 0, het$Q_pval[1], NA),
      egger_intercept = ifelse(!is.null(pleio) && nrow(pleio) > 0, pleio$egger_intercept[1], NA),
      egger_p = ifelse(!is.null(pleio) && nrow(pleio) > 0, pleio$pval[1], NA),
      presso_global_p = presso_global,
      presso_outlier = presso_outlier,
      note = "OK",
      stringsAsFactors = FALSE
    )

    # 追加其他方法
    for (m in c("MR Egger", "Weighted median", "Weighted mode")) {
      mrow <- mr_res[mr_res$method == m, ]
      if (nrow(mrow) > 0) {
        res_df <- rbind(res_df, data.frame(
          exposure = exp$name, outcome = outc$name,
          exposure_label = exp$label, outcome_label = outc$label,
          n_snps = n_iv, method = m,
          b = mrow$b, se = mrow$se, pval = mrow$pval,
          or = exp(mrow$b), or_lci = exp(mrow$b - 1.96 * mrow$se),
          or_uci = exp(mrow$b + 1.96 * mrow$se),
          fstat = round(fstat, 1), steiger = steiger,
          het_q = NA, het_p = NA, egger_intercept = NA, egger_p = NA,
          presso_global_p = NA, presso_outlier = NA,
          note = "OK",
          stringsAsFactors = FALSE
        ))
      }
    }

    # 日志
    cat(sprintf("  IVW: OR=%.3f (%.3f-%.3f) P=%.2e | F=%.1f | Steiger=%s | nSNP=%d\n",
                exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), pval,
                fstat, steiger, n_iv),
        file = log_file, append = TRUE)

    res_df
  }
}

# ---- 错误/跳过时的空行 ----
df_error <- function(exp, outc, n_iv = NA, note = "error") {
  data.frame(
    exposure = exp$name, outcome = outc$name,
    exposure_label = exp$label, outcome_label = outc$label,
    n_snps = n_iv, method = "IVW",
    b = NA, se = NA, pval = NA,
    or = NA, or_lci = NA, or_uci = NA,
    fstat = NA, steiger = NA,
    het_q = NA, het_p = NA,
    egger_intercept = NA, egger_p = NA,
    presso_global_p = NA, presso_outlier = NA,
    note = note,
    stringsAsFactors = FALSE
  )
}

# ---- IVs 缓存文件路径 ----
iv_cache_file <- "data/iv_cache.rds"
results_csv <- "results/tables/03_batch_mr_results.csv"
ivs_csv <- "results/tables/03_batch_ivs.csv"

# ---- 加载已有结果（支持断点续跑） ----
if (file.exists(results_csv)) {
  existing_results <- read.csv(results_csv, stringsAsFactors = FALSE)
  done_pairs <- unique(paste(existing_results$exposure, existing_results$outcome, sep = "->"))
  cat("已有结果:", length(done_pairs), "对，跳过\n\n", file = log_file, append = TRUE)
} else {
  existing_results <- NULL
  done_pairs <- character(0)
}

# IV 缓存
iv_cache <- list()
if (file.exists(iv_cache_file)) {
  iv_cache <- readRDS(iv_cache_file)
}

# ---- 主循环 ----
total_pairs <- length(exposures) * length(outcomes)
pair_num <- 0
all_ivs <- list()

for (exp in exposures) {
  # 提取该暴露的 IVs（缓存）
  if (exp$id %in% names(iv_cache)) {
    exp_dat <- iv_cache[[exp$id]]
    cat("IV缓存命中:", exp$label, "(", nrow(exp_dat), "IVs )\n", file = log_file, append = TRUE)
  } else {
    cat("提取IVs:", exp$label, "...\n", file = log_file, append = TRUE)
    exp_dat <- safe_extract_instruments(exp$id, p1 = 5e-8)

    if (inherits(exp_dat, "error") || is.null(exp_dat) || nrow(exp_dat) < 3) {
      cat("  P<5e-8 不够，试 P<1e-5...\n", file = log_file, append = TRUE)
      exp_dat <- safe_extract_instruments(exp$id, p1 = 1e-5)
    }

    if (inherits(exp_dat, "error") || is.null(exp_dat) || nrow(exp_dat) < 2) {
      cat("  SKIP:", exp$label, "IVs不足\n", file = log_file, append = TRUE)
      # 为该暴露的所有结局写错误行
      for (outc in outcomes) {
        pair_num <- pair_num + 1
        pair_id <- paste(exp$name, outc$name, sep = "->")
        if (pair_id %in% done_pairs) next
        err_row <- df_error(exp, outc, n_iv = 0, note = "insufficient_IVs")
        if (is.null(existing_results)) {
          existing_results <- err_row
        } else {
          existing_results <- rbind(existing_results, err_row)
        }
        write.csv(existing_results, results_csv, row.names = FALSE)
      }
      next
    }

    iv_cache[[exp$id]] <- exp_dat
    saveRDS(iv_cache, iv_cache_file)
    cat("  IVs:", nrow(exp_dat), "(已缓存)\n", file = log_file, append = TRUE)
  }

  # 遍历结局
  for (outc in outcomes) {
    pair_num <- pair_num + 1
    pair_id <- paste(exp$name, outc$name, sep = "->")

    # 断点续跑：跳过已完成的
    if (pair_id %in% done_pairs) {
      cat(sprintf("[%d/%d] %s -> %s (跳过，已完成)\n", pair_num, total_pairs,
                  exp$label, outc$label),
          file = log_file, append = TRUE)
      next
    }

    msg <- sprintf("[%d/%d] %s (%s) -> %s (%s)",
                   pair_num, total_pairs, exp$label, exp$id, outc$label, outc$id)
    cat(msg, "\n", file = log_file, append = TRUE)
    cat(msg, "\n")

    # 运行分析（tryCatch 在调用层，不在函数内用 return）
    result <- tryCatch(
      run_mr_pair(exp, outc, exp_dat),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n", file = log_file, append = TRUE)
        df_error(exp, outc, n_iv = nrow(exp_dat),
                 note = paste("ERROR:", substr(conditionMessage(e), 1, 50)))
      }
    )

    # 收集 IVs
    if (!is.null(result) && result$note[1] == "OK") {
      # 重新提取 harmonised dat 用于 IV 记录
      out_dat <- safe_extract_outcome(exp_dat$SNP, outc$id)
      if (!inherits(out_dat, "error") && !is.null(out_dat)) {
        dat <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
        if (!is.null(dat)) {
          dat <- subset(dat, mr_keep == TRUE)
          if (nrow(dat) > 0) {
            all_ivs[[pair_id]] <- data.frame(
              exposure = exp$name, outcome = outc$name,
              SNP = dat$SNP, beta.exposure = dat$beta.exposure,
              se.exposure = dat$se.exposure, pval.exposure = dat$pval.exposure,
              beta.outcome = dat$beta.outcome, se.outcome = dat$se.outcome,
              pval.outcome = dat$pval.outcome,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    # 逐对追加保存（不丢进度）
    if (is.null(existing_results)) {
      write.csv(result, results_csv, row.names = FALSE)
      existing_results <- result
    } else {
      existing_results <- rbind(existing_results, result)
      write.csv(existing_results, results_csv, row.names = FALSE)
    }

    cat("  已保存 ->", results_csv, "\n", file = log_file, append = TRUE)
  }
}

# ---- 保存 IVs ----
if (length(all_ivs) > 0) {
  all_ivs_df <- do.call(rbind, all_ivs)
  write.csv(all_ivs_df, ivs_csv, row.names = FALSE)
}

# ---- IVW 汇总 ----
ivw_only <- existing_results[existing_results$method == "IVW" |
                               existing_results$method == "Inverse variance weighted", ]
cat("\n\n=== IVW 结果汇总 ===\n", file = log_file, append = TRUE)
for (i in seq_len(nrow(ivw_only))) {
  r <- ivw_only[i, ]
  sig <- ifelse(!is.na(r$pval) && r$pval < 0.05, "*", "")
  cat(sprintf("%s -> %s: OR=%.3f (%.3f-%.3f) P=%.3e %s [nSNP=%s F=%s Steiger=%s] %s\n",
              r$exposure_label, r$outcome_label,
              r$or, r$or_lci, r$or_uci, r$pval, sig,
              as.character(r$n_snps), as.character(r$fstat), as.character(r$steiger),
              r$note),
      file = log_file, append = TRUE)
}

cat("\n=== 批量 MR v2 完成 ===\n", file = log_file, append = TRUE)
cat("结果:", results_csv, "\n")
cat("IVs:", ivs_csv, "\n")
cat("=== 完成 ===\n")
