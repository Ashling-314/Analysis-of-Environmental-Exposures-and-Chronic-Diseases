# Leave-one-out for 14 IVW-significant pairs (with correct ID mapping)
suppressMessages(suppressWarnings({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
}))

jwt_token <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(jwt_token)) {
  rl <- readLines("C:/Users/[username]/Documents/.Renviron", warn = FALSE)
  jl <- grep("^OPENGWAS_JWT=", rl, value = TRUE)
  jwt_token <- sub("^OPENGWAS_JWT=", "", jl[1])
}
Sys.setenv(OPENGWAS_JWT = jwt_token)
options(timeout = 120)

# ID mappings
exp_id_map <- c(
  "Smoking_initiation" = "ieu-b-4877",
  "Alcohol_intake_freq" = "ukb-b-5779",
  "NO2_air_pollution" = "ukb-b-9942",
  "Coffee_intake" = "ukb-b-5237",
  "Fresh_fruit_intake" = "ukb-b-10817"
)

# Outcome ID map (from the original batch MR script)
out_id_map <- c(
  "CAD" = "ieu-a-7",
  "LungCancer" = "ieu-a-966",
  "IschStroke" = "ebi-a-GCST005843",
  "HeartFailure" = "ebi-a-GCST009541",
  "COPD" = "ebi-a-GCST90018807",
  "T2D" = "ebi-a-GCST006867"
)

# Load IV cache
iv_cache <- readRDS("data/iv_cache.rds")

# Load significant pairs
p2 <- fread("results/tables/03_batch_mr_results.csv")
p2[, pval := as.numeric(pval)]
sig <- p2[method == "IVW" & pval < 0.05]

cat("Will process", nrow(sig), "pairs for LOO\n\n")

loo_summary <- data.table()

for (i in 1:nrow(sig)) {
  row <- sig[i]
  exp_name <- row$exposure
  out_name <- row$outcome

  exp_id <- exp_id_map[exp_name]
  out_id <- out_id_map[out_name]

  if (is.na(exp_id) || is.na(out_id)) {
    cat(sprintf("[%d] SKIP: %s->%s (no ID mapping)\n", i, exp_name, out_name))
    next
  }

  cat(sprintf("[%d/%d] %s -> %s (%s -> %s)\n", i, nrow(sig),
              row$exposure_label, row$outcome_label, exp_id, out_id))

  exposure_dat <- iv_cache[[exp_id]]
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
    cat("  No cached exposure IVs\n")
    next
  }

  tryCatch({
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = out_id,
      proxies = FALSE
    )
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      cat("  No outcome data\n")
      next
    }

    dat <- harmonise_data(exposure_dat, outcome_dat)
    dat <- as.data.frame(dat)
    dat <- dat[!is.na(dat$beta.outcome) & !is.na(dat$beta.exposure), ]
    if (nrow(dat) < 3) {
      cat("  Too few SNPs:", nrow(dat), "\n")
      next
    }

    # Leave-one-out
    loo <- mr_leaveoneout(dat, method = mr_ivw)
    loo_dt <- as.data.table(loo)

    # Check stability: does removing any SNP change direction or make p > 0.05?
    main_beta <- row$b
    main_p <- row$pval
    flipped <- sum(sign(loo_dt$b) != sign(main_beta))
    n_loo_p_gt_05 <- sum(loo_dt$p > 0.05)

    loo_summary <- rbind(loo_summary, data.table(
      pair = paste(row$exposure_label, "->", row$outcome_label),
      n_snps = nrow(dat),
      main_beta = main_beta,
      main_p = main_p,
      loo_min_beta = min(loo_dt$b),
      loo_max_beta = max(loo_dt$b),
      loo_min_p = min(loo_dt$p),
      loo_max_p = max(loo_dt$p),
      n_flipped = flipped,
      n_loo_p_gt_05 = n_loo_p_gt_05,
      stable = (flipped == 0 & n_loo_p_gt_05 == 0)
    ))

    cat(sprintf("  LOO done: stable=%s, flipped=%d, n_p>0.05=%d\n",
                flipped == 0 & n_loo_p_gt_05 == 0, flipped, n_loo_p_gt_05))
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })

  Sys.sleep(2)
}

fwrite(loo_summary, "results/tables/08_leaveoneout_summary.csv")

cat("\n\n=== LOO 稳定性总结 ===\n")
cat("总对数:", nrow(loo_summary), "\n")
cat("完全稳定 (无方向翻转, 无p>0.05):", sum(loo_summary$stable), "/", nrow(loo_summary), "\n")
cat("有方向翻转:", sum(loo_summary$n_flipped > 0), "/", nrow(loo_summary), "\n")
cat("至少一个LOO p>0.05:", sum(loo_summary$n_loo_p_gt_05 > 0), "/", nrow(loo_summary), "\n\n")

if (nrow(loo_summary[!stable]) > 0) {
  cat("--- 不稳定的对 (需注意) ---\n")
  for (i in 1:nrow(loo_summary[!stable])) {
    r <- loo_summary[!stable][i]
    cat(sprintf("  %s: 翻转=%d, p>0.05=%d\n", r$pair, r$n_flipped, r$n_loo_p_gt_05))
  }
}

cat("\n结果已保存: results/tables/08_leaveoneout_summary.csv\n")
