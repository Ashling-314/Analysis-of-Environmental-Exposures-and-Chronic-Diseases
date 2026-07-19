# ============================================================
# 07_coloc.R — Phase 7 共定位分析 (v4 健壮版)
# 改进: 超时控制 + 断点续跑 + 暴露数据复用 + 限速规避
# ============================================================

suppressMessages(suppressWarnings({
  library(ieugwasr)
  library(coloc)
  library(data.table)
}))

# ---- 关键: 设置 API 调用超时, 防止无限卡死 ----
options(timeout = 120)          # 单次 HTTP 请求最多 120 秒
# 不设置 ieugwasr_api, 使用默认 URL (https://gwas-api.mrcieu.ac.uk/)

# ---- 设置 JWT 令牌 ----
jwt_token <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(jwt_token)) {
  renviron_file <- "C:/Users/WuYan/Documents/.Renviron"
  if (file.exists(renviron_file)) {
    rl <- readLines(renviron_file, warn = FALSE)
    jl <- grep("^OPENGWAS_JWT=", rl, value = TRUE)
    if (length(jl) > 0) jwt_token <- sub("^OPENGWAS_JWT=", "", jl[1])
  }
}
Sys.setenv(OPENGWAS_JWT = jwt_token)
cat("JWT length:", nchar(jwt_token), "\n")

# ---- 日志函数 (带 flush) ----
log_file <- "logs/07_coloc.log"
log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_file, append = TRUE)
  cat(msg, "\n")
  flush.console()
}

# ---- 断点续跑: 读取已有结果, 跳过已完成的对 ----
done_pairs <- character(0)
results_csv <- "results/tables/07_coloc_results.csv"
if (file.exists(results_csv)) {
  old <- tryCatch(fread(results_csv), error = function(e) NULL)
  if (!is.null(old) && nrow(old) > 0) {
    done_pairs <- paste(old$exposure, "->", old$outcome)
    log("已存在结果 ", nrow(old), " 对, 将跳过: ", paste(done_pairs, collapse = "; "))
  }
}

log("=== Phase 7 共定位开始 (v4) ===")
log("时间:", format(Sys.time()))
log("JWT length:", nchar(jwt_token))

# ---- 加载 1000G .bim ----
bim_file <- "D:/kindle_ideas/kindle_ideas/idea15_线粒体DNA异质性/data/ref/1000G_EUR/EUR.bim"
log("加载 1000G .bim...")
bim <- fread(bim_file, header = FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
bim[, chr := as.integer(chr)]
bim[, pos := as.numeric(pos)]
bim[, cm := NULL]
bim[, a1 := NULL]
bim[, a2 := NULL]
gc()
log("  .bim 总 SNP 数:", nrow(bim))

# ---- 获取 lead SNP (硬编码已知 lead SNP, 避免 tophits() 间歇性 HTML 错误) ----
hardcoded_leads <- list(
  "ieu-b-4877" = list(chr = 11, pos = 112911004, snp = "rs7938812", p = 1e-300),
  "ukb-b-5237" = list(chr = 15, pos = 75027880,  snp = "rs2472297", p = 1e-50),
  "ukb-b-5779" = list(chr = 4,  pos = 100239319, snp = "rs1229984", p = 1e-100)
)

get_lead_snp <- function(id) {
  # 优先使用硬编码的 lead SNP
  if (id %in% names(hardcoded_leads)) {
    lead <- hardcoded_leads[[id]]
    log("  使用硬编码 lead SNP: ", lead$snp, " (chr", lead$chr, ":", lead$pos, ")")
    return(lead)
  }
  # 回退到 tophits()
  th <- tryCatch(tophits(id = id), error = function(e) {
    log("  tophits error: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(th) || nrow(th) == 0) return(NULL)
  th <- th[order(th$p), ]
  chr_val <- as.integer(th$chr[1])
  pos_val <- as.numeric(th$position[1])
  if (is.na(chr_val) || is.na(pos_val)) {
    log("  WARNING: chr 或 pos 为 NA!")
    return(NULL)
  }
  return(list(chr = chr_val, pos = pos_val, snp = th$rsid[1], p = th$p[1]))
}

# ---- 获取区域统计量 (带超时+重试) ----
get_regional_stats <- function(id, target_chr, target_pos, window = 1.5e5) {
  chr_match <- bim$chr == target_chr
  pos_lo <- target_pos - window
  pos_hi <- target_pos + window
  region_mask <- chr_match & (bim$pos >= pos_lo & bim$pos <= pos_hi)

  n_region <- sum(region_mask, na.rm = TRUE)
  log("    区域: chr", target_chr, ":", pos_lo, "-", pos_hi,
      " | SNP 数:", n_region)

  if (n_region == 0) {
    log("    无区域 SNP, SKIP")
    return(NULL)
  }

  region_snps <- bim$rsid[region_mask]

  # 分批查询, 每批 300 SNP, 批间 sleep 3s
  all_res <- list()
  batch_size <- 300
  n_batches <- ceiling(n_region / batch_size)
  log("    分批查询: ", n_batches, " 批 (每批 ", batch_size, " SNP)")

  for (i in seq(1, n_region, by = batch_size)) {
    batch_idx <- i:min(i + batch_size - 1, n_region)
    batch <- region_snps[batch_idx]
    batch_num <- ceiling(i / batch_size)

    res <- NULL
    # 最多重试 2 次
    for (attempt in 1:3) {
      res <- tryCatch({
        associations(variants = batch, id = id)
      }, error = function(e) {
        err_msg <- conditionMessage(e)
        log("    batch ", batch_num, "/", n_batches,
            " 尝试", attempt, " error: ", substr(err_msg, 1, 80))
        NULL
      })
      if (!is.null(res) && nrow(res) > 0) break
      if (attempt < 3) Sys.sleep(5 * attempt)  # 递增等待
    }

    if (!is.null(res) && nrow(res) > 0) {
      all_res[[length(all_res) + 1]] <- res
    }
    if (batch_num %% 3 == 0 || batch_num == n_batches) {
      log("    batch ", batch_num, "/", n_batches, " done (成功 ",
          length(all_res), "/", batch_num, ")")
    }
    Sys.sleep(3)  # 批间睡眠, 避免 API 限速
  }

  if (length(all_res) == 0) {
    log("    所有批次失败, SKIP")
    return(NULL)
  }
  result <- do.call(rbind, all_res)
  log("    总获取 SNP 数:", nrow(result))
  return(result)
}

# ---- 构造 coloc 数据框 ----
make_coloc_df <- function(reg, default_N) {
  snp_col <- intersect(c("rsid", "SNP", "snp"), names(reg))[1]
  beta_col <- intersect(c("beta", "b"), names(reg))[1]
  se_col <- intersect(c("se", "SE"), names(reg))[1]
  p_col <- intersect(c("p", "pval"), names(reg))[1]
  eaf_col <- intersect(c("eaf", "MAF", "maf"), names(reg))
  eaf_col <- if (length(eaf_col) > 0) eaf_col[1] else NULL
  n_col <- intersect(c("n", "N", "samplesize"), names(reg))
  n_col <- if (length(n_col) > 0) n_col[1] else NULL

  df <- data.frame(
    snp = reg[[snp_col]],
    beta = as.numeric(reg[[beta_col]]),
    varbeta = as.numeric(reg[[se_col]])^2,
    pvalues = as.numeric(reg[[p_col]]),
    N = if (!is.null(n_col)) as.numeric(reg[[n_col]]) else default_N,
    MAF = if (!is.null(eaf_col)) as.numeric(reg[[eaf_col]]) else NA,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$snp) & !is.na(df$beta) & !is.na(df$varbeta) & df$varbeta > 0, ]
  if (all(is.na(df$MAF))) {
    log("    WARNING: 无 MAF 数据, 将跳过该区域")
    return(NULL)
  }
  df$MAF <- ifelse(df$MAF > 0.5, 1 - df$MAF, df$MAF)
  df <- df[!is.na(df$MAF) & df$MAF > 0, ]
  return(df)
}

# ---- 共定位对定义 ----
coloc_pairs <- list(
  list(exp_id = "ieu-b-4877",  exp_name = "Smoking",  exp_type = "cc",
       out_id = "ebi-a-GCST90018807", out_name = "COPD",      out_type = "cc",
       exp_N = 607291, out_N = 4616),
  list(exp_id = "ieu-b-4877",  exp_name = "Smoking",  exp_type = "cc",
       out_id = "ieu-a-966",         out_name = "LungCancer", out_type = "cc",
       exp_N = 607291, out_N = 11269),
  list(exp_id = "ieu-b-4877",  exp_name = "Smoking",  exp_type = "cc",
       out_id = "ieu-a-7",           out_name = "CAD",        out_type = "cc",
       exp_N = 607291, out_N = 86995),
  list(exp_id = "ukb-b-5237",  exp_name = "Coffee",   exp_type = "quant",
       out_id = "ebi-a-GCST006867",  out_name = "T2D",        out_type = "cc",
       exp_N = 913819, out_N = 65566),
  list(exp_id = "ukb-b-5779",  exp_name = "Alcohol",  exp_type = "quant",
       out_id = "ebi-a-GCST90018807", out_name = "COPD",      out_type = "cc",
       exp_N = 913819, out_N = 4616)
)

# ---- 主循环 ----
all_coloc <- list()

# 载入已有结果到 all_coloc (用于断点续跑时的最终汇总)
if (file.exists(results_csv)) {
  old <- tryCatch(fread(results_csv), error = function(e) NULL)
  if (!is.null(old) && nrow(old) > 0) {
    all_coloc <- lapply(seq_len(nrow(old)), function(i) old[i, ])
  }
}

log("\n开始共定位分析...")
log("已完成对: ", ifelse(length(done_pairs) > 0, paste(done_pairs, collapse = "; "), "无"))

# 暴露数据缓存 (同一 exp_id 只查一次)
exp_cache <- list()

for (pi in seq_along(coloc_pairs)) {
  pair <- coloc_pairs[[pi]]
  pair_key <- paste(pair$exp_name, "->", pair$out_name)

  if (pair_key %in% done_pairs) {
    log("\n[跳过] Pair ", pi, "/", length(coloc_pairs), ": ", pair_key, " (已完成)")
    next
  }

  log("\n========================================")
  log("Pair ", pi, "/", length(coloc_pairs), ": ", pair_key)
  log("========================================")

  # 获取 lead SNP
  lead <- get_lead_snp(pair$exp_id)
  if (is.null(lead)) {
    log("  无 lead SNP, SKIP")
    next
  }
  log("  Lead SNP:", lead$snp, " chr", lead$chr, " pos", lead$pos,
      " p=", format(lead$p, digits = 3))

  chr_val <- lead$chr
  pos_val <- lead$pos
  lead_snp <- lead$snp

  # 获取暴露区域数据 (带缓存)
  exp_reg <- NULL
  if (pair$exp_id %in% names(exp_cache)) {
    log("  [1] 复用暴露区域数据缓存 (", pair$exp_id, ")")
    exp_reg <- exp_cache[[pair$exp_id]]
  } else {
    log("  [1] 获取暴露区域数据 (", pair$exp_id, ")")
    exp_reg <- get_regional_stats(pair$exp_id, chr_val, pos_val)
    if (!is.null(exp_reg) && nrow(exp_reg) > 0) {
      exp_cache[[pair$exp_id]] <- exp_reg
    }
  }
  if (is.null(exp_reg) || nrow(exp_reg) == 0) {
    log("  无暴露区域数据, SKIP")
    next
  }

  # 获取结局区域数据
  log("  [2] 获取结局区域数据 (", pair$out_id, ")")
  out_reg <- get_regional_stats(pair$out_id, chr_val, pos_val)
  if (is.null(out_reg) || nrow(out_reg) == 0) {
    log("  无结局区域数据, SKIP")
    next
  }

  # 构造 coloc 数据框
  log("  [3] 构造 coloc 数据框")
  exp_df <- make_coloc_df(exp_reg, pair$exp_N)
  out_df <- make_coloc_df(out_reg, pair$out_N)

  if (is.null(exp_df) || is.null(out_df)) {
    log("  coloc 数据框构造失败, SKIP")
    next
  }

  common <- intersect(exp_df$snp, out_df$snp)
  log("  交集 SNP 数:", length(common),
      " (exp:", nrow(exp_df), ", out:", nrow(out_df), ")")

  if (length(common) < 50) {
    log("  交集 SNP < 50, SKIP (区域数据不足)")
    next
  }

  exp_df <- exp_df[exp_df$snp %in% common, ]
  out_df <- out_df[out_df$snp %in% common, ]
  exp_df <- exp_df[order(exp_df$snp), ]
  out_df <- out_df[out_df$snp %in% common, ]
  out_df <- out_df[order(out_df$snp), ]

  # 运行 coloc
  log("  [4] 运行 coloc.abf")
  D1 <- list(
    snp = exp_df$snp, beta = exp_df$beta, varbeta = exp_df$varbeta,
    pvalues = exp_df$pvalues, N = exp_df$N, MAF = exp_df$MAF,
    type = pair$exp_type
  )
  D2 <- list(
    snp = out_df$snp, beta = out_df$beta, varbeta = out_df$varbeta,
    pvalues = out_df$pvalues, N = out_df$N, MAF = out_df$MAF,
    type = pair$out_type
  )

  res <- tryCatch(coloc.abf(dataset1 = D1, dataset2 = D2), error = function(e) {
    log("  coloc error: ", conditionMessage(e))
    NULL
  })

  if (is.null(res)) next

  pp <- res$summary
  log("  PP.H4=", round(pp["PP.H4.abf"], 4),
      " PP.H3=", round(pp["PP.H3.abf"], 4),
      " n_snps=", length(common))

  result_row <- data.frame(
    exposure = pair$exp_name, outcome = pair$out_name,
    chr = chr_val, pos = pos_val, lead_snp = lead_snp,
    n_snps = length(common),
    PP_H0 = pp["PP.H0.abf"], PP_H1 = pp["PP.H1.abf"],
    PP_H2 = pp["PP.H2.abf"], PP_H3 = pp["PP.H3.abf"],
    PP_H4 = pp["PP.H4.abf"], coloc_5percent = pp["PP.H4.abf"] >= 0.5,
    note = "OK", stringsAsFactors = FALSE
  )

  all_coloc[[length(all_coloc) + 1]] <- result_row
  write.csv(do.call(rbind, all_coloc), results_csv, row.names = FALSE)
  log("  结果已保存 (累计:", length(all_coloc), ")")

  # pair 间 sleep
  if (pi < length(coloc_pairs)) Sys.sleep(5)
}

log("\n=== Phase 7 共定位完成 ===")
log("时间:", format(Sys.time()))
log("总共定位区域数:", length(all_coloc))
if (length(all_coloc) > 0) {
  final <- do.call(rbind, all_coloc)
  log("PP.H4 >= 0.5 的区域:", sum(final$PP_H4 >= 0.5, na.rm = TRUE))
  write.csv(final, results_csv, row.names = FALSE)
}
