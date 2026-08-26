#' Score and classify clusters against CeNGEN cell types
#'
#' Convenience wrapper that runs [score_cengen_matrix()] followed by
#' [classify_cengen_matches()], and attaches the full cluster x cell-type
#' score matrix to the result (retrieve it with [cengen_matrix()]) for
#' drilling into all candidate cell types beyond the top 2.
#'
#' @inheritParams score_cengen_matrix
#' @param cutoff,min_gap,min_genes_used Passed to [classify_cengen_matches()].
#'
#' @return A `cengen_scores` tibble, one row per cluster (see
#'   [classify_cengen_matches()]), with the full score matrix attached as
#'   an attribute.
#' @export
score_cengen_clusters <- function(
  markers,
  reference,
  top_n = 20,
  min_pct_expr = 10,
  combine = c("geometric", "arithmetic"),
  cov_weight = 0.5,
  spec_weight = 0.5,
  cutoff = 0.6,
  min_gap = 0.10,
  min_genes_used = 5,
  cluster_col = "cluster",
  gene_col = "gene"
) {
  combine <- match.arg(combine)

  matrix <- score_cengen_matrix(
    markers, reference,
    top_n = top_n, min_pct_expr = min_pct_expr,
    combine = combine, cov_weight = cov_weight, spec_weight = spec_weight,
    cluster_col = cluster_col, gene_col = gene_col
  )

  result <- classify_cengen_matches(
    matrix,
    cutoff = cutoff, min_gap = min_gap, min_genes_used = min_genes_used
  )

  attr(result, "matrix") <- matrix
  result
}

#' Get the full cluster x cell-type score matrix
#'
#' Retrieves the score matrix attached to a `cengen_scores` object by
#' [score_cengen_clusters()].
#'
#' @param x A `cengen_scores` object.
#' @return The attached long score matrix (see [score_cengen_matrix()]), or
#'   `NULL` if none is attached.
#' @export
cengen_matrix <- function(x) {
  attr(x, "matrix")
}
