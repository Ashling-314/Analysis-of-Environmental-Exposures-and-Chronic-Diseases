# Coloc with narrower ±50kb window
suppressMessages(suppressWarnings({
  library(ieugwasr)
  library(coloc)
  library(data.table)
}))

jwt_token <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(jwt_token)) {
  rl <- readLines("C:/Users/WuYan/Documents/.Renviron", warn = FALSE)
  jl <- grep("^OPENGWAS_JWT=", rl, value = TRUE)
  jwt_token <- sub("^OPENGWAS_JWT=", "", jl[1])
}
Sys.setenv(OPENGWAS_JWT = jwt_token)
options(timeout = 120)

log_file <- "logs/09_coloc_narrow.log"
log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_file, append = TRUE)
  cat(msg, "\n")
  flush.console()
}

log("=== 共定位 ±50kb 窗口重跑 ===")
log("时间:", format(Sys.time()))

# Same hardcoded lead SNPs
hardcoded_leads <- list(
  "ieu-b-4877" = list(chr = 11, pos = 112911004, snp = "rs7938812"),
  "ukb-b-5237" = list(chr = 15, pos = 75027880,  snp = "rs2472297"),
  "ukb-b-5779" = list(chr = 4,  pos = 100239319, snp = "rs1229984")
)

# Pairs to run (same 3 remaining)
pairs_to_run <- list(
  list(exp_id = "ukb-b-5237", exp_name = "Coffee",   exp_type = "quant",
       out_id = "ebi-a-GCST006867",  out_name = "T2D",        out_type = "cc",
       exp_N = 913819, out_N = 65566),
  list(exp_id = "ukb-b-5779", exp_name = "Alcohol",  exp_type = "quant",
       out_id = "ebi-a-GCST90018807", out_name = "COPD",      out_type = "cc",
       exp_N = 913819, out_N = 468475),
  list(exp_id = "ieu-b-4877", exp_name = "Smoking",  exp_type = "cc",
       out_id = "ieu-a-7",           out_name = "CAD",        out_type = "cc",
       exp_N = 607291, out_N = 86995)
)

# Load pre-filtered .bim (only 3 regions, ~1934 SNPs)
log("加载预过滤 .bim (3个区域, 1934 SNP)...")
bim <- fread("data/1000G_3regions.bim", header = FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
bim[, chr := as.integer(chr)]
bim[, pos := as.numeric(pos)]
bim[, cm := NULL][, a1 := NULL][, a2 := NULL]
gc()
log("  SNP 数:", nrow(bim))

# Get lead SNP
get_lead <- function(id) {
  if (id %in% names(hardcoded_leads)) return(hardcoded_leads[[id]])
  return(NULL)
}

# Get regional stats (same logic as v4 but with 5e4 window)
get_regional <- function(id, target_chr, target_pos, window = 5e4) {
  region_mask <- bim$chr == target_chr & bim$pos >= (target_pos - window) & bim$pos <= (target_pos + window)
  region_snps <- bim$rsid[region_mask]
  n_region <- length(region_snps)
  log("    区域: chr", target_chr, ":", target_pos - window, "-", target_pos + window,
      " | SNPs:", n_region)
  if (n_region == 0) return(NULL)

  all_res <- list()
  batch_size <- 300
  n_batches <- ceiling(n_region / batch_size)
  for (batch_num in 1:n_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, n_region)
    batch_snps <- region_snps[start_idx:end_idx]
    for (attempt in 1:3) {
      res <- tryCatch(associations(variants = batch_snps, id = id), error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        all_res[[batch_num]] <- res
        break
      }
      Sys.sleep(5)
    }
    Sys.sleep(3)
  }
  do.call(rbind, all_res)
}

# Coloc wrapper
make_coloc_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  col_map <- c("p" = "pval", "beta" = "beta", "se" = "se", "MAF" = "eaf",
               "n" = "n", "pos" = "pos", "snp" = "rsid", "chr" = "chr")
  for (old in names(col_map)) {
    if (old %in% colnames(df) && !col_map[old] %in% colnames(df)) {
      df[[col_map[old]]] <- df[[old]]
    }
  }
  if (!"MAF" %in% colnames(df) && "eaf" %in% colnames(df)) {
    df$MAF <- pmin(df$eaf, 1 - df$eaf)
  }
  if (!"N" %in% colnames(df) && "n" %in% colnames(df)) {
    df$N <- df$n
  }
  if (!"varbeta" %in% colnames(df) && "se" %in% colnames(df)) {
    df$varbeta <- df$se^2
  }
  return(df)
}

# Output
results_csv <- "results/tables/09_coloc_narrow_50kb.csv"
if (!file.exists(results_csv)) {
  cat("exposure,outcome,chr,pos,lead_snp,n_snps,PP_H0,PP_H1,PP_H2,PP_H3,PP_H4,coloc_5percent,note\n",
      file = results_csv)
}

for (pair in pairs_to_run) {
  log("\n========================================")
  log("Pair: ", pair$exp_name, " -> ", pair$out_name)
  lead <- get_lead(pair$exp_id)
  if (is.null(lead)) { log("  No lead SNP, SKIP"); next }
  log("  Lead: ", lead$snp, " chr", lead$chr, ":", lead$pos)

  # Exposure
  log("  [1] Exposure regional stats...")
  exp_res <- get_regional(pair$exp_id, lead$chr, lead$pos)
  if (is.null(exp_res)) { log("  No exposure data"); next }
  exp_df <- make_coloc_df(exp_res)
  exp_df <- exp_df[!is.na(exp_df$beta) & !is.na(exp_df$varbeta) & exp_df$varbeta > 0, ]
  log("    Got", nrow(exp_df), "SNPs")

  # Outcome
  log("  [2] Outcome regional stats...")
  out_res <- get_regional(pair$out_id, lead$chr, lead$pos)
  if (is.null(out_res)) { log("  No outcome data"); next }
  out_df <- make_coloc_df(out_res)
  out_df <- out_df[!is.na(out_df$beta) & !is.na(out_df$varbeta) & out_df$varbeta > 0, ]
  log("    Got", nrow(out_df), "SNPs")

  # Common SNPs
  common <- intersect(exp_df$rsid, out_df$rsid)
  log("  Common SNPs:", length(common))
  if (length(common) < 10) { log("  Too few common SNPs"); next }

  exp_sub <- exp_df[exp_df$rsid %in% common, ]
  out_sub <- out_df[out_df$rsid %in% common, ]
  exp_sub <- exp_sub[match(common, exp_sub$rsid), ]
  out_sub <- out_sub[match(common, out_sub$rsid), ]

  exp_N <- pair$exp_N
  out_N <- pair$out_N

  # Coloc
  log("  [3] Running coloc.abf...")
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      coloc.abf(
        dataset1 = list(beta = exp_sub$beta, varbeta = exp_sub$varbeta,
                        MAF = exp_sub$MAF, type = pair$exp_type, N = exp_N, snp = exp_sub$rsid),
        dataset2 = list(beta = out_sub$beta, varbeta = out_sub$varbeta,
                        MAF = out_sub$MAF, type = pair$out_type, N = out_N, snp = out_sub$rsid)
      )
    )),
    error = function(e) { log("  coloc error:", conditionMessage(e)); NULL }
  )
  if (is.null(res)) next

  pp <- res$summary
  log("  PP.H4 =", round(pp["PP.H4.abf"], 4))

  line <- sprintf('"%s","%s",%d,%d,"%s",%d,%.6e,%.6e,%.6e,%.6e,%.6e,%s,"%s"',
                  pair$exp_name, pair$out_name, lead$chr, lead$pos, lead$snp, length(common),
                  pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"],
                  ifelse(pp["PP.H4.abf"] >= 0.5, "TRUE", "FALSE"), "OK_50kb")
  cat(line, "\n", file = results_csv, append = TRUE)
}

log("\n=== 50kb 共定位完成 ===")
log("时间:", format(Sys.time()))
