#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
workspace_root <- "."

required <- c("highmean", "highDmean", "matrixStats")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

processed_dir <- file.path(workspace_root, "data")
index_dir <- file.path(workspace_root, "data", "indices")
result_dir <- Sys.getenv("LEUKEMIA_RESULT_DIR", file.path(workspace_root, "results", "eight_methods"))
metadata_dir <- file.path(workspace_root, "results", "metadata")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

read_matrix_gz <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  dat <- tryCatch(utils::read.csv(con, check.names = FALSE), finally = close(con))
  features <- setdiff(names(dat), c("sample_id", "group"))
  x <- as.matrix(dat[, features, drop = FALSE])
  storage.mode(x) <- "double"
  list(data = x, sample_id = as.character(dat$sample_id), group = as.character(dat$group), feature_names = features)
}

mixed <- read_matrix_gz(file.path(processed_dir, "leukemia_mixed250.csv.gz"))
top <- read_matrix_gz(file.path(processed_dir, "leukemia_top250.csv.gz"))

hdbf_file <- file.path(workspace_root, "R", "HDBF_simulation_study.R")
if (!file.exists(hdbf_file)) stop("Missing frozen HDBF reference: ", hdbf_file)
hdbf <- new.env(parent = globalenv())
sys.source(hdbf_file, envir = hdbf)

run_eight <- function(x, y) {
  clx <- highmean::apval_Cai2014(x, y, eq.cov = FALSE)
  cq <- highmean::apval_Chen2010(x, y, eq.cov = FALSE)
  cctl <- hdbf$cauchy_combination(clx$pval, cq$pval)
  zwl <- highDmean::zwl_test(x, y, order = 0)
  skk <- highDmean::SKK_test(x, y)
  wcct <- hdbf$wcct_pvalue(x, y)
  h1 <- tcrossprod(x)
  h2 <- tcrossprod(y)
  tn <- hdbf$chen_qin_Tn(x, y, h1, h2)
  variance <- max(
    1e-12,
    (2 / (nrow(x) * (nrow(x) - 1))) * hdbf$trace_square_est(h1, nrow(x)) +
      (2 / (nrow(y) * (nrow(y) - 1))) * hdbf$trace_square_est(h2, nrow(y)) +
      (4 / (nrow(x) * nrow(y))) * hdbf$cross_trace_est(x, y)
  )
  qn <- tn / sqrt(variance)
  ewcct <- hdbf$ewcct_pvalue(wcct, qn, nrow(x), nrow(y), ncol(x), threshold_coef = 2.1)
  ewccts <- hdbf$cauchy_combination(ewcct, skk$pvalue)
  data.frame(
    CLX = as.numeric(clx$pval), CQ = as.numeric(cq$pval), CCTL = as.numeric(cctl),
    ZWLm = as.numeric(zwl$pvalue), SKK = as.numeric(skk$pvalue), WCCT = as.numeric(wcct),
    EWCCT = as.numeric(ewcct), EWCCTS = as.numeric(ewccts),
    Qn = as.numeric(qn), SKK_statistic = as.numeric(skk$TSvalue),
    stringsAsFactors = FALSE
  )
}

top_all <- which(top$group == "ALL")
top_aml <- which(top$group == "AML")
real <- run_eight(top$data[top_all, , drop = FALSE], top$data[top_aml, , drop = FALSE])
real_long <- data.frame(
  method = c("CLX", "CQ", "CCTL", "ZWLm", "SKK", "WCCT", "EWCCT", "EWCCTS"),
  p_value = as.numeric(real[1, c("CLX", "CQ", "CCTL", "ZWLm", "SKK", "WCCT", "EWCCT", "EWCCTS")]),
  stringsAsFactors = FALSE
)
utils::write.csv(real_long, file.path(result_dir, "realdata_top250.csv"), row.names = FALSE, quote = FALSE)
utils::write.csv(real, file.path(result_dir, "realdata_top250_auxiliary.csv"), row.names = FALSE, quote = FALSE)

null_path <- file.path(index_dir, "leukemia_null_indices_1based.csv")
alt_path <- file.path(index_dir, "leukemia_alternative_indices_1based.csv")
if (!file.exists(null_path) || !file.exists(alt_path)) stop("Run R/02_reproduce_cq.R first to create shared indices.")
null_indices <- utils::read.csv(null_path, check.names = FALSE)
alt_indices <- utils::read.csv(alt_path, check.names = FALSE)
use_reps <- Sys.getenv("LEUKEMIA_REPS", "")
if (nzchar(use_reps)) {
  use_reps <- as.integer(use_reps)
  if (is.na(use_reps) || use_reps < 1L) stop("LEUKEMIA_REPS must be a positive integer when set.")
  null_indices <- head(null_indices, use_reps)
  alt_indices <- head(alt_indices, use_reps)
}

extract_pair <- function(dat, row, scenario) {
  if (scenario == "null_h0") {
    x_idx <- as.integer(row[grep("^g1_", names(row))])
    y_idx <- as.integer(row[grep("^g2_", names(row))])
  } else {
    x_idx <- as.integer(row[grep("^all_", names(row))])
    y_idx <- as.integer(row[grep("^aml_", names(row))])
  }
  list(x = dat$data[x_idx, , drop = FALSE], y = dat$data[y_idx, , drop = FALSE])
}

run_replicates <- function(dat, indices, scenario) {
  out <- vector("list", nrow(indices))
  for (b in seq_len(nrow(indices))) {
    pair <- extract_pair(dat, indices[b, , drop = FALSE], scenario)
    z <- run_eight(pair$x, pair$y)
    out[[b]] <- data.frame(replicate = b, scenario = scenario, z, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

replicates <- rbind(
  run_replicates(mixed, null_indices, "null_h0"),
  run_replicates(mixed, alt_indices, "alternative")
)
utils::write.csv(replicates, file.path(result_dir, "resampling_replicates_r.csv"), row.names = FALSE, quote = FALSE)

methods <- c("CLX", "CQ", "CCTL", "ZWLm", "SKK", "WCCT", "EWCCT", "EWCCTS")
summary_rows <- vector("list", length(methods) * 2L)
counter <- 0L
for (scenario in c("null_h0", "alternative")) {
  z <- replicates[replicates$scenario == scenario, , drop = FALSE]
  for (method in methods) {
    counter <- counter + 1L
    p <- z[[method]]
    valid <- is.finite(p)
    reject <- valid & p < 0.05
    count <- sum(reject)
    total <- sum(valid)
    rate <- if (total) count / total else NA_real_
    ci <- if (total) stats::binom.test(count, total)$conf.int else c(NA_real_, NA_real_)
    summary_rows[[counter]] <- data.frame(
      scenario = scenario, method = method, alpha = 0.05,
      repetitions = length(p), valid_repetitions = total, rejects = count,
      rate = rate, monte_carlo_se = if (total) sqrt(rate * (1 - rate) / total) else NA_real_,
      ci95_lower = ci[1], ci95_upper = ci[2], stringsAsFactors = FALSE
    )
  }
}
utils::write.csv(do.call(rbind, summary_rows), file.path(result_dir, "resampling_summary_r.csv"), row.names = FALSE, quote = FALSE)

writeLines(c(
  paste0("reference_md5=", unname(tools::md5sum(hdbf_file))),
  "covariance_for_CLX_CQ=unequal (eq.cov=FALSE), matching HDBF Behrens-Fisher workflow",
  "ZWLm=highDmean::zwl_test(order=0)",
  "SKK=highDmean::SKK_test package pvalue, two-sided as used by HDBF_simulation_study.R",
  "EWCCTS=EWCCT combined with the same two-sided SKK package pvalue",
  "threshold_coef=2.1"
), file.path(metadata_dir, "eight_methods_r.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(metadata_dir, "sessionInfo_eight_methods_r.txt"))

cat("Eight-method real-data output written to ", result_dir, ".\n", sep = "")
