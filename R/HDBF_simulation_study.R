# ==============================================================================
# Simulation study for the HDBF manuscript
#
# Manuscript:
#   Power-Enhanced Cauchy Combination Tests for the Two-Sample
#   High-Dimensional Behrens-Fisher Problem
#
# Purpose:
#   Reproduce the Monte Carlo simulation workflow for type-I-error and power
#   summaries of CLX, CQ, CCTL, ZWLm, SKK, WCCT, EWCCT, and EWCCTS
#   under AR(1) and block covariance structures.
#
# Tested with:
#   R 4.4.3 on Windows
#
# Main entry points:
#   source("HDBF_simulation_study.R")
#   results <- run_all()
#   type1 <- run_type1_error()
#   power <- run_power()
#
# Notes:
#   - The full simulation may take substantial time because it uses 2000 Monte
#     Carlo replications per covariance structure by default.
#   - To run a smaller smoke check, modify hdbf_config() after sourcing this
#     file, for example cfg <- hdbf_config(); cfg$n_simulations <- 10.
# ==============================================================================

#### Package requirements ####

required_packages <- c(
  "future","future.apply",
  "highmean","highDmean",
  "matrixStats","mvnfast")

# Stop early if a package required by the simulation is not installed.
check_packages <- function(packages = required_packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0L) {
    stop(
      "Please install the following R packages before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Return the default simulation settings used in the submitted experiment.
hdbf_config <- function() {
  list(
    n_simulations = 2000L,
    alpha = 0.05,
    n1 = 30L,
    n2 = 30L,
    d = 300L,
    structures = c("AR1", "Block"),
    rho1 = 0.3,
    rho2 = 0.3,
    sparsity = 0.01,
    power_signal_C = 3.0,
    threshold_coefs = c(2.1),
    num_blocks = 20L,
    workers = max(1L, as.integer(parallel::detectCores()) - 4L, na.rm = TRUE),
    use_progress = TRUE,
    output_dir = getwd()
  )
}

# Run an expression under the requested future parallel plan, then restore it.
with_parallel_plan <- function(config, expr) {
  check_packages()

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (config$workers > 1L) {
    future::plan(future::multisession, workers = config$workers)
  } else {
    future::plan(future::sequential)
  }

  force(expr)
}

#### Test statistics and p-value combinations ####

# Compute the Chen-Qin quadratic-form statistic using precomputed Gram matrices.
chen_qin_Tn <- function(X1, X2, H1, H2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  term1 <- (sum(H1) - sum(diag(H1))) / (n1 * (n1 - 1))
  term2 <- (sum(H2) - sum(diag(H2))) / (n2 * (n2 - 1))
  term3 <- -2 * sum(tcrossprod(X1, X2)) / (n1 * n2)
  term1 + term2 + term3
}

# Estimate tr(Sigma^2) from a Gram matrix for one sample.
trace_square_est <- function(H, m) {
  if (m < 4L) {
    return(NA_real_)
  }

  H_sq_elem <- H^2
  row_sums_H <- rowSums(H)
  row_sums_H_sq <- rowSums(H_sq_elem)
  HH <- H %*% H

  S_sq_all <- outer(row_sums_H_sq, row_sums_H_sq, "+") - 2 * HH
  S_val_all <- outer(row_sums_H, row_sums_H, "-")

  diag_H <- diag(H)
  Di <- matrix(diag_H, m, m)
  Dj <- t(Di)

  diff_at_i <- Di - H
  diff_at_j <- H - Dj

  S_sq_valid <- S_sq_all - (diff_at_i^2) - (diff_at_j^2)
  S_val_valid <- S_val_all - diff_at_i - diff_at_j

  n_valid <- m - 2L
  inner_sum <- 2 * (n_valid * S_sq_valid - S_val_valid^2)
  diag(inner_sum) <- 0

  sum(inner_sum) / (4 * m * (m - 1) * (m - 2) * (m - 3))
}

# Estimate the cross trace tr(Sigma_1 Sigma_2) from two centered samples.
cross_trace_est <- function(X, Y) {
  m <- nrow(X)
  n <- nrow(Y)
  Xc <- t(t(X) - matrixStats::colMeans2(X))
  Yc <- t(t(Y) - matrixStats::colMeans2(Y))
  cross_terms <- tcrossprod(Yc, Xc)
  sum(cross_terms^2) / ((m - 1) * (n - 1))
}

# Compute the WCCT p-value from coordinate-wise Welch t-tests.
wcct_pvalue <- function(X1, X2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  m1 <- matrixStats::colMeans2(X1)
  m2 <- matrixStats::colMeans2(X2)
  v1 <- matrixStats::colVars(X1)
  v2 <- matrixStats::colVars(X2)

  vn1 <- v1 / n1
  vn2 <- v2 / n2
  se <- sqrt(vn1 + vn2)
  se[se < 1e-12] <- 1e-12

  t_stat <- (m1 - m2) / se
  df_num <- (vn1 + vn2)^2
  df_den <- (vn1^2) / (n1 - 1) + (vn2^2) / (n2 - 1)
  df <- df_num / df_den

  p_vals <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  tan_vals <- tan((0.5 - p_vals) * pi)
  tan_vals[!is.finite(tan_vals) | tan_vals > 1e10] <- 1e10

  0.5 - atan(mean(tan_vals)) / pi
}

# Combine two p-values using the Cauchy combination rule.
cauchy_combination <- function(p1, p2) {
  if (is.na(p1) || is.na(p2)) {
    return(NA_real_)
  }
  t1 <- tan((0.5 - p1) * pi)
  t2 <- tan((0.5 - p2) * pi)
  if (!is.finite(t1)) {
    t1 <- 1e10
  }
  if (!is.finite(t2)) {
    t2 <- 1e10
  }
  pcauchy(0.5 * t1 + 0.5 * t2, lower.tail = FALSE)
}

# Compute the EWCCT p-values using the finite enhancement term in the manuscript.
ewcct_pvalue <- function(p_wcct, Qn, n1, n2, d, threshold_coef = 2.1) {
  p_ew <- p_wcct
  triggered <- !is.na(p_wcct) & !is.na(Qn) & abs(Qn) > sqrt(threshold_coef * log(d))

  if (any(triggered)) {
    T_wcct <- tan((0.5 - p_wcct[triggered]) * pi)
    T_ewcct <- T_wcct + n1 + n2
    p_ew[triggered] <- 0.5 - atan(T_ewcct) / pi
  }

  p_ew
}

#### Data generation ####

# Build the covariance matrix for one simulation scenario.
covariance_matrix <- function(d, structure, rho, num_blocks = NULL) {
  if (structure == "AR1") {
    return(rho^abs(outer(seq_len(d), seq_len(d), "-")))
  }

  if (structure == "Block") {
    if (is.null(num_blocks) || is.na(num_blocks)) {
      stop("num_blocks is required for the Block covariance structure.", call. = FALSE)
    }
    block_size <- ceiling(d / num_blocks)
    Sigma <- matrix(0, nrow = d, ncol = d)

    for (block_id in seq_len(num_blocks)) {
      start_idx <- (block_id - 1L) * block_size + 1L
      end_idx <- min(block_id * block_size, d)
      if (start_idx > d) {
        break
      }
      idx <- start_idx:end_idx
      Sigma[idx, idx] <- rho
      diag(Sigma)[idx] <- 1
    }

    diag(Sigma) <- diag(Sigma) + 1e-6
    return(Sigma)
  }

  stop("Unknown covariance structure: ", structure, call. = FALSE)
}

# Define the sparse mean shift used in the power experiment.
mean_shift <- function(config, signal_C) {
  if (signal_C == 0) {
    return(list(k = 0L, value = 0))
  }

  k <- max(1L, ceiling(config$d * config$sparsity))
  effective_n <- config$n1 * config$n2 / (config$n1 + config$n2)
  value <- signal_C * sqrt(log(config$d) / effective_n) / sqrt(k)
  list(k = k, value = value)
}

# Generate the two samples for one Monte Carlo replication.
simulate_samples <- function(config, chol1, chol2, signal_C) {
  signal <- mean_shift(config, signal_C)

  mu1 <- rep(0, config$d)
  mu2 <- rep(0, config$d)
  if (signal$k > 0L) {
    mu2[sample(seq_len(config$d), signal$k)] <- signal$value
  }

  X1 <- mvnfast::rmvn(config$n1, mu1, chol1, isChol = TRUE)
  X2 <- mvnfast::rmvn(config$n2, mu2, chol2, isChol = TRUE)

  list(X1 = X1, X2 = X2, k = signal$k, value = signal$value)
}

#### Simulation ####

# Run one Monte Carlo replication and return the statistics used for summaries.
one_replication <- function(config, chol1, chol2, equal_covariance, signal_C) {
  samples <- simulate_samples(config, chol1, chol2, signal_C)
  X1 <- samples$X1
  X2 <- samples$X2

  H1 <- tcrossprod(X1)
  H2 <- tcrossprod(X2)

  p_wcct <- wcct_pvalue(X1, X2)

  Tn <- chen_qin_Tn(X1, X2, H1, H2)
  tr_s1 <- trace_square_est(H1, config$n1)
  tr_s2 <- trace_square_est(H2, config$n2)
  tr_s12 <- cross_trace_est(X1, X2)

  var_est <- max(
    1e-12,
    (2 / (config$n1 * (config$n1 - 1))) * tr_s1 +
      (2 / (config$n2 * (config$n2 - 1))) * tr_s2 +
      (4 / (config$n1 * config$n2)) * tr_s12
  )
  Qn <- Tn / sqrt(var_est)

  p_cq <- tryCatch(
    highmean::apval_Chen2010(X1, X2, eq.cov = equal_covariance)$pval,
    error = function(e) NA_real_
  )
  p_clx <- tryCatch(
    highmean::apval_Cai2014(X1, X2, eq.cov = equal_covariance)$pval,
    error = function(e) NA_real_
  )
  p_zwl <- tryCatch(
    highDmean::zwl_test(X1, X2)$pvalue,
    error = function(e) NA_real_
  )
  p_skk <- tryCatch(
    highDmean::SKK_test(X1, X2)$pvalue,
    error = function(e) NA_real_
  )

  c(
    wcct = p_wcct,
    Qn = Qn,
    clx = p_clx,
    cq = p_cq,
    zwl = p_zwl,
    skk = p_skk
  )
}

# Run all replications for one scenario, with optional console progress.
run_replications <- function(config, chol1, chol2, equal_covariance, signal_C) {
  indices <- seq_len(config$n_simulations)

  if (!isTRUE(config$use_progress)) {
    return(future.apply::future_lapply(
      indices,
      function(idx) {
        one_replication(config, chol1, chol2, equal_covariance, signal_C)
      },
      future.packages = required_packages
    ))
  }

  chunk_size <- max(1L, min(length(indices), config$workers * 5L))
  chunks <- split(indices, ceiling(seq_along(indices) / chunk_size))
  sim_results <- vector("list", length(indices))
  progress_bar <- utils::txtProgressBar(
    min = 0,
    max = length(indices),
    initial = 0,
    style = 3
  )
  on.exit(close(progress_bar), add = TRUE)

  completed <- 0L
  for (chunk in chunks) {
    chunk_results <- future.apply::future_lapply(
      chunk,
      function(idx) {
        one_replication(config, chol1, chol2, equal_covariance, signal_C)
      },
      future.seed = TRUE,
      future.packages = required_packages
    )
    sim_results[chunk] <- chunk_results
    completed <- completed + length(chunk)
    utils::setTxtProgressBar(progress_bar, completed)
  }
  cat("\n")
  sim_results
}

# Prepare covariance matrices and run the replications for one structure.
run_scenario <- function(structure, config, signal_C) {
  S1 <- covariance_matrix(config$d, structure, config$rho1, config$num_blocks)
  S2 <- covariance_matrix(config$d, structure, config$rho2, config$num_blocks)
  chol1 <- chol(S1)
  chol2 <- chol(S2)
  equal_covariance <- isTRUE(all.equal(S1, S2, tolerance = 1e-10, check.attributes = FALSE))

  if (isTRUE(config$use_progress)) {
    cat(sprintf(
      "Running %s condition (signal C = %s, simulations = %d)\n",
      structure,
      format(signal_C, trim = TRUE, scientific = FALSE),
      config$n_simulations
    ))
  }

  sim_results <- run_replications(config, chol1, chol2, equal_covariance, signal_C)
  res_mat <- as.data.frame(do.call(rbind, sim_results))
  res_mat[] <- lapply(res_mat, as.numeric)
  res_mat
}

#### Summaries and output ####

# Estimate the rejection rate at the configured significance level.
rejection_rate <- function(x, alpha) {
  valid <- !is.na(x)
  if (!any(valid)) {
    return(NA_real_)
  }
  mean(x[valid] < alpha)
}

# Convert replication-level p-values into method-level rejection rates.
summarise_rates <- function(res_mat, config, prefix) {
  rates <- c(
    CLX = rejection_rate(res_mat$clx, config$alpha),
    CQ = rejection_rate(res_mat$cq, config$alpha),
    CCTL = rejection_rate(mapply(cauchy_combination, res_mat$clx, res_mat$cq), config$alpha),
    ZWLm = rejection_rate(res_mat$zwl, config$alpha),
    SKK = rejection_rate(res_mat$skk, config$alpha),
    WCCT = rejection_rate(res_mat$wcct, config$alpha)
  )

  for (coef in config$threshold_coefs) {
    p_ew <- ewcct_pvalue(res_mat$wcct, res_mat$Qn, config$n1, config$n2, config$d, coef)
    rates["EWCCT"] <- rejection_rate(p_ew, config$alpha)
    rates["EWCCTS"] <- rejection_rate(
      mapply(cauchy_combination, p_ew, res_mat$skk),
      config$alpha
    )
  }

  names(rates) <- paste(prefix, names(rates), sep = "_")
  as.list(rates)
}

# Summarise one covariance structure as one row of the results table.
scenario_summary <- function(structure, config, signal_C, prefix) {
  scenario_results <- run_scenario(
    structure = structure,
    config = config,
    signal_C = signal_C
  )

  metadata <- list(
    structure = structure,
    n1 = config$n1,
    n2 = config$n2,
    d = config$d,
    rho1 = config$rho1,
    rho2 = config$rho2,
    alpha = config$alpha,
    n_simulations = config$n_simulations
  )

  if (prefix == "power") {
    metadata$sparsity <- config$sparsity
  }

  c(metadata, summarise_rates(scenario_results, config, prefix))
}

# Bind a list of named rows into a plain data frame.
rows_to_data_frame <- function(rows) {
  df <- as.data.frame(do.call(rbind, lapply(rows, as.data.frame)), stringsAsFactors = FALSE)
  numeric_cols <- setdiff(names(df), "structure")
  df[numeric_cols] <- lapply(df[numeric_cols], function(x) as.numeric(as.character(x)))
  rownames(df) <- NULL
  df
}

# Write a result table to CSV and store the output path as an attribute.
write_results <- function(results, config, result_type) {
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
  filename <- sprintf(
    "HDBF_%s_results_%s.csv",
    result_type,
    format(Sys.time(), "%Y%m%d_%H%M")
  )
  output_file <- file.path(config$output_dir, filename)
  write.csv(results, output_file, row.names = FALSE)
  attr(results, "output_file") <- output_file
  cat(sprintf("Results saved to: %s\n", output_file))
  results
}

#### Public entry points ####

# Run the type-I-error simulation under the null mean shift.
run_type1_error <- function(config = hdbf_config(), write_csv = TRUE) {
  with_parallel_plan(config, {
    rows <- lapply(
      config$structures,
      function(structure) {
        scenario_summary(
          structure,
          config,
          signal_C = 0,
          prefix = "type1"
        )
      }
    )
    results <- rows_to_data_frame(rows)
    if (isTRUE(write_csv)) {
      results <- write_results(results, config, "type1_error")
    }
    results
  })
}

# Run the power simulation using the configured sparse mean shift.
run_power <- function(config = hdbf_config(), write_csv = TRUE) {
  with_parallel_plan(config, {
    rows <- lapply(
      config$structures,
      function(structure) {
        scenario_summary(
          structure,
          config,
          signal_C = config$power_signal_C,
          prefix = "power"
        )
      }
    )
    results <- rows_to_data_frame(rows)
    if (isTRUE(write_csv)) {
      results <- write_results(results, config, "power")
    }
    results
  })
}

# Run both type-I-error and power simulations.
run_all <- function(config = hdbf_config(), write_csv = TRUE) {
  type1_error <- run_type1_error(config, write_csv = write_csv)
  power <- run_power(config, write_csv = write_csv)
  results <- list(type1_error = type1_error, power = power)

  if (isTRUE(write_csv)) {
    attr(results, "output_files") <- c(
      type1_error = attr(type1_error, "output_file"),
      power = attr(power, "output_file")
    )
  }
  results
}

#### Main ####
# Run the full submission simulation when this file is executed directly.
main <- function() {
  config <- hdbf_config()
  cat("HDBF simulation\n")
  cat(sprintf("n1 = %d, n2 = %d, d = %d\n", config$n1, config$n2, config$d))
  cat(sprintf("structures = %s\n", paste(config$structures, collapse = ", ")))
  cat(sprintf("rho1 = %.3f, rho2 = %.3f\n", config$rho1, config$rho2))
  cat(sprintf("n_simulations = %d, workers = %d\n", config$n_simulations, config$workers))
  start_time <- proc.time()
  run_all(config, write_csv = TRUE)
  elapsed <- (proc.time() - start_time)[3]
  cat(sprintf("Completed in %.2f seconds.\n", elapsed))
}

if (sys.nframe() == 0L) {
  main()
}
