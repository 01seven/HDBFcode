#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Run from the workspace root.  Using getwd() avoids Windows code-page
# conversion of the Chinese workspace path when Rscript receives --file.
workspace_root <- "."

required <- c("golubEsets", "Biobase")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

processed_dir <- file.path(workspace_root, "data")
metadata_dir <- file.path(workspace_root, "results", "metadata")
validation_dir <- file.path(workspace_root, "results", "validation")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, quote = FALSE)
}

stable_order <- function(p) {
  # Original probe order is the deterministic tie breaker.
  order(p, seq_along(p), na.last = TRUE, method = "radix")
}

safe_t_pvalue <- function(x, y, var_equal) {
  out <- tryCatch(
    stats::t.test(x, y, var.equal = var_equal, paired = FALSE)$p.value,
    error = function(e) NA_real_
  )
  as.numeric(out)
}

data(Golub_Merge, package = "golubEsets")
expr <- Biobase::exprs(Golub_Merge)
pheno <- Biobase::pData(Golub_Merge)

if (!identical(dim(expr), c(7129L, 72L))) {
  stop("Golub_Merge must be 7,129 genes x 72 samples; got ", paste(dim(expr), collapse = " x "))
}
if (!all(is.finite(expr))) stop("Golub_Merge contains NA/Inf values.")
if (!"ALL.AML" %in% colnames(pheno)) stop("Golub_Merge lacks pData column ALL.AML.")
group <- as.character(pheno$ALL.AML)
if (!identical(as.integer(table(group)[c("ALL", "AML")]), c(47L, 25L))) {
  stop("Expected ALL=47 and AML=25; got ", paste(names(table(group)), table(group), collapse = "; "))
}

feature_id <- rownames(expr)
if (is.null(feature_id) || anyDuplicated(feature_id)) stop("Feature IDs must be present and unique.")
sample_id <- colnames(expr)
if (is.null(sample_id) || anyDuplicated(sample_id)) stop("Sample IDs must be present and unique.")

all_idx <- which(group == "ALL")
aml_idx <- which(group == "AML")
welch_p <- vapply(seq_len(nrow(expr)), function(i) {
  safe_t_pvalue(expr[i, all_idx], expr[i, aml_idx], var_equal = FALSE)
}, numeric(1))
pooled_p <- vapply(seq_len(nrow(expr)), function(i) {
  safe_t_pvalue(expr[i, all_idx], expr[i, aml_idx], var_equal = TRUE)
}, numeric(1))

welch_rank <- stable_order(welch_p)
pooled_rank <- stable_order(pooled_p)
top250_idx <- welch_rank[seq_len(250L)]
mixed_idx <- c(welch_rank[seq_len(50L)], welch_rank[(length(welch_rank) - 199L):length(welch_rank)])
if (length(unique(mixed_idx)) != 250L) stop("mixed250 selection is not unique.")

ranking_welch <- data.frame(
  rank = seq_along(welch_rank),
  original_index = welch_rank,
  feature_id = feature_id[welch_rank],
  p_value = welch_p[welch_rank],
  selected_top250 = welch_rank %in% top250_idx,
  selected_mixed250 = welch_rank %in% mixed_idx
)
ranking_pooled <- data.frame(
  rank = seq_along(pooled_rank),
  original_index = pooled_rank,
  feature_id = feature_id[pooled_rank],
  p_value = pooled_p[pooled_rank],
  selected_top250 = pooled_rank %in% pooled_rank[seq_len(250L)],
  selected_mixed250 = pooled_rank %in% c(pooled_rank[seq_len(50L)], pooled_rank[(length(pooled_rank) - 199L):length(pooled_rank)])
)
utils::write.csv(ranking_welch, file.path(processed_dir, "gene_ranking_welch.csv"), row.names = FALSE, quote = FALSE)
utils::write.csv(ranking_pooled, file.path(processed_dir, "gene_ranking_pooled.csv"), row.names = FALSE, quote = FALSE)

sample_metadata <- data.frame(
  sample_position = seq_along(sample_id),
  sample_id = sample_id,
  group = group,
  ALL_or_AML_index = NA_integer_,
  stringsAsFactors = FALSE
)
sample_metadata$ALL_or_AML_index <- ave(seq_along(group), group, FUN = seq_along)
utils::write.csv(sample_metadata, file.path(processed_dir, "sample_metadata.csv"), row.names = FALSE, quote = FALSE)

# Preserve the exact package expression values and the two fixed 250-feature matrices.
expr_samples <- as.data.frame(t(expr), check.names = FALSE)
colnames(expr_samples) <- feature_id
expr_samples <- cbind(sample_id = sample_id, group = group, expr_samples)
write_csv_gz(expr_samples, file.path(processed_dir, "leukemia_expression_72x7129.csv.gz"))

make_matrix_frame <- function(indices) {
  z <- as.data.frame(t(expr[indices, , drop = FALSE]), check.names = FALSE)
  colnames(z) <- feature_id[indices]
  cbind(sample_id = sample_id, group = group, z)
}
write_csv_gz(make_matrix_frame(top250_idx), file.path(processed_dir, "leukemia_top250.csv.gz"))
write_csv_gz(make_matrix_frame(mixed_idx), file.path(processed_dir, "leukemia_mixed250.csv.gz"))

hash_paths <- c(
  expression_72x7129 = file.path(processed_dir, "leukemia_expression_72x7129.csv.gz"),
  top250 = file.path(processed_dir, "leukemia_top250.csv.gz"),
  mixed250 = file.path(processed_dir, "leukemia_mixed250.csv.gz"),
  ranking_welch = file.path(processed_dir, "gene_ranking_welch.csv"),
  ranking_pooled = file.path(processed_dir, "gene_ranking_pooled.csv"),
  sample_metadata = file.path(processed_dir, "sample_metadata.csv")
)
writeLines(paste(names(hash_paths), unname(tools::md5sum(hash_paths))), file.path(metadata_dir, "processed_md5.txt"), useBytes = TRUE)

overlap_top250 <- length(intersect(welch_rank[seq_len(250L)], pooled_rank[seq_len(250L)]))
overlap_mixed <- length(intersect(
  c(welch_rank[seq_len(50L)], welch_rank[(length(welch_rank) - 199L):length(welch_rank)]),
  c(pooled_rank[seq_len(50L)], pooled_rank[(length(pooled_rank) - 199L):length(pooled_rank)])
))
preprocessing_audit <- data.frame(
  item = c("source", "dimensions", "ALL", "AML", "extra_transform", "range_min", "range_max", "missing", "welch_p_top250", "pooled_p_top250", "welch_p_mixed_top50_cut", "welch_p_mixed_bottom200_cut", "ranking_overlap_top250", "ranking_overlap_mixed250"),
  value = c(
    "golubEsets::Golub_Merge", paste(dim(expr), collapse = "x"), sum(group == "ALL"), sum(group == "AML"),
    "none; direct package exprs", min(expr), max(expr), sum(!is.finite(expr)),
    min(welch_p[top250_idx]), min(pooled_p[pooled_rank[seq_len(250L)]]),
    max(welch_p[welch_rank[seq_len(50L)]]), min(welch_p[welch_rank[(length(welch_rank) - 199L):length(welch_rank)]]),
    overlap_top250, overlap_mixed
  ), stringsAsFactors = FALSE
)
utils::write.csv(preprocessing_audit, file.path(validation_dir, "preprocessing_audit.csv"), row.names = FALSE, quote = FALSE)

provenance <- c(
  paste0("source=golubEsets::Golub_Merge"),
  paste0("golubEsets_version=", as.character(utils::packageVersion("golubEsets"))),
  paste0("Biobase_version=", as.character(utils::packageVersion("Biobase"))),
  paste0("dimensions=", paste(dim(expr), collapse = "x")),
  paste0("groups=ALL:", sum(group == "ALL"), ";AML:", sum(group == "AML")),
  "ranking_primary=two-sided Welch t-test, stable original-probe-order tie break",
  "ranking_sensitivity=two-sided pooled-variance t-test",
  paste0("top250_features=", paste(feature_id[top250_idx], collapse = ",")),
  paste0("mixed250_features=", paste(feature_id[mixed_idx], collapse = ",")),
  paste0("R=", getRversion())
)
writeLines(provenance, file.path(metadata_dir, "provenance_prepare.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(metadata_dir, "sessionInfo_prepare.txt"))

cat("Prepared Golub_Merge: ", paste(dim(expr), collapse = " x "), "\n", sep = "")
cat("ALL/AML: ", sum(group == "ALL"), "/", sum(group == "AML"), "\n", sep = "")
cat("Welch top250 minimum p: ", format(min(welch_p[top250_idx]), digits = 8), "\n", sep = "")
cat("Welch/pooled top250 overlap: ", overlap_top250, "/250\n", sep = "")
