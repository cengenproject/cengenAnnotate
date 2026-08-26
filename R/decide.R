#' Classify clusters as annotatable or in need of review
#'
#' Applies the cutoff/gap decision rule to a cluster x cell-type score
#' matrix (from [score_cengen_matrix()]), reducing each cluster to its best
#' and second-best matching cell type and a decision.
#'
#' Decision priority, evaluated in order:
#' \enumerate{
#'   \item If `n_genes_used < min_genes_used`: `"review_insufficient_markers"`,
#'     regardless of score (too few usable marker genes to trust any match).
#'   \item Else if `best_score >= cutoff` AND `gap >= min_gap`:
#'     `"annotate"` (one cell type stands out clearly).
#'   \item Otherwise: `"review_ambiguous"` (either no candidate scores high
#'     enough, or more than one candidate is plausible).
#' }
#'
#' @param score_matrix A long score matrix, as returned by
#'   [score_cengen_matrix()].
#' @param cutoff Minimum `best_score` required to auto-annotate.
#' @param min_gap Minimum gap between the best and second-best score
#'   required to auto-annotate (guards against ambiguous ties).
#' @param min_genes_used Minimum number of usable marker genes required;
#'   below this, the cluster is always sent to review.
#'
#' @return A tibble with one row per cluster: `cluster`, `best_cell_type`,
#'   `best_score`, `second_cell_type`, `second_score`, `gap`,
#'   `n_genes_used`, `n_genes_requested`, `n_genes_missing_from_reference`,
#'   `n_genes_uninformative`, `decision`. Sorted by `best_score` descending
#'   (most confident cluster first), not by input cluster order.
#' @export
classify_cengen_matches <- function(
  score_matrix,
  cutoff = 0.6,
  min_gap = 0.10,
  min_genes_used = 5
) {
  assert_has_columns(
    score_matrix,
    c(
      "cluster", "cell_type", "score", "n_genes_used", "n_genes_requested",
      "n_genes_missing_from_reference", "n_genes_uninformative"
    ),
    "score_matrix", "cengen_invalid_score_matrix_error"
  )

  clusters <- unique(score_matrix$cluster)

  rows <- lapply(clusters, function(cl) {
    sub <- score_matrix[score_matrix$cluster == cl, , drop = FALSE]
    ord <- order(-ifelse(is.na(sub$score), -Inf, sub$score), sub$cell_type)
    sub <- sub[ord, , drop = FALSE]

    best_cell_type <- sub$cell_type[1]
    best_score <- sub$score[1]
    second_cell_type <- if (nrow(sub) >= 2) sub$cell_type[2] else NA_character_
    second_score <- if (nrow(sub) >= 2) sub$score[2] else NA_real_
    gap <- best_score - second_score

    n_genes_used <- sub$n_genes_used[1]

    eps <- 1e-9
    decision <- if (is.na(n_genes_used) || n_genes_used < min_genes_used) {
      "review_insufficient_markers"
    } else if (!is.na(best_score) && !is.na(gap) &&
      best_score >= cutoff - eps && gap >= min_gap - eps) {
      "annotate"
    } else {
      "review_ambiguous"
    }

    tibble::tibble(
      cluster = cl,
      best_cell_type = best_cell_type,
      best_score = best_score,
      second_cell_type = second_cell_type,
      second_score = second_score,
      gap = gap,
      n_genes_used = n_genes_used,
      n_genes_requested = sub$n_genes_requested[1],
      n_genes_missing_from_reference = sub$n_genes_missing_from_reference[1],
      n_genes_uninformative = sub$n_genes_uninformative[1],
      decision = decision
    )
  })

  result <- dplyr::bind_rows(rows)
  # Sort by confidence (best_score, descending) rather than input cluster
  # order, so the clusters worth trusting most are at the top and the
  # murkiest calls sink to the bottom. NA scores (all-excluded marker
  # panels) sort last.
  result <- result[order(-ifelse(is.na(result$best_score), -Inf, result$best_score)), , drop = FALSE]
  class(result) <- c("cengen_scores", class(result))
  result
}
