#' Combine cengenAnnotate results with Seurat TransferData predictions
#'
#' Augments a `cengen_scores` object (from [score_cengen_clusters()]) with
#' per-cluster predictions from `Seurat::TransferData()`, producing a
#' 3-tier classification validated across two independent datasets:
#' clusters `cengenAnnotate` already confidently called stay `"annotate"`;
#' clusters it didn't, but where TransferData's own per-cell votes are
#' confident, get promoted to `"provisional"` (labeled with TransferData's
#' call); everything else is left as-is for manual review.
#'
#' On both test datasets, `cengenAnnotate`'s `"annotate"` calls and
#' TransferData's confident per-cell-vote majority never disagreed. A
#' genuine conflict (both methods confident, different answers) is
#' therefore treated as a rare edge case rather than silently resolved:
#' it's flagged as `"review_conflict"`, keeping `cengenAnnotate`'s own
#' (independently validated) call as `best_cell_type` while recording
#' TransferData's call in `td_call` for inspection.
#'
#' `object` must already have `Seurat::TransferData()` results attached,
#' e.g.:
#' ```r
#' anchors <- Seurat::FindTransferAnchors(reference = ref, query = object)
#' predictions <- Seurat::TransferData(anchorset = anchors, refdata = ref$cell_type)
#' object <- Seurat::AddMetaData(object, metadata = predictions)
#' ```
#'
#' @param scores A `cengen_scores` object, as returned by
#'   [score_cengen_clusters()].
#' @param object A Seurat object with TransferData results already added
#'   (i.e. has `predicted_col` and `score_col` metadata columns).
#' @param cluster_col Name of the cluster identity column in
#'   `object@meta.data` (default `"seurat_clusters"`).
#' @param predicted_col,score_col Names of the per-cell TransferData
#'   columns (Seurat's own defaults: `"predicted.id"`,
#'   `"prediction.score.max"`).
#' @param td_min_purity Minimum fraction of a cluster's cells that must
#'   agree on TransferData's majority vote for that cluster to count as
#'   TransferData-confident.
#' @param td_min_score Minimum mean per-cell `prediction.score.max` within
#'   a cluster for it to count as TransferData-confident.
#'
#' @return A `cengen_scores` tibble (same shape as `scores`, with the score
#'   matrix attribute preserved so [cengen_matrix()] still works), with
#'   `decision` extended to include `"provisional"` and
#'   `"review_conflict"`, `best_cell_type` updated to TransferData's call
#'   for `"provisional"` rows, and three new diagnostic columns: `td_call`,
#'   `td_purity`, `td_mean_score`.
#' @export
combine_with_transfer_data <- function(
  scores,
  object,
  cluster_col = "seurat_clusters",
  predicted_col = "predicted.id",
  score_col = "prediction.score.max",
  td_min_purity = 0.6,
  td_min_score = 0.6
) {
  if (!inherits(object, "Seurat")) {
    cengen_abort("`object` must be a Seurat object.", class = "cengen_invalid_object_error")
  }
  assert_has_columns(
    object@meta.data, c(cluster_col, predicted_col, score_col),
    "object@meta.data", "cengen_invalid_object_error"
  )
  assert_has_columns(
    scores, c("cluster", "best_cell_type", "best_score", "decision"),
    "scores", "cengen_invalid_scores_error"
  )

  td <- tibble::tibble(
    cluster = as.character(object@meta.data[[cluster_col]]),
    predicted = as.character(object@meta.data[[predicted_col]]),
    pred_score = as.numeric(object@meta.data[[score_col]])
  )

  td_summary <- dplyr::summarise(
    dplyr::group_by(td, .data$cluster),
    td_call = names(sort(table(.data$predicted), decreasing = TRUE))[1],
    td_purity = max(table(.data$predicted)) / dplyr::n(),
    td_mean_score = mean(.data$pred_score, na.rm = TRUE),
    .groups = "drop"
  )

  matrix_attr <- attr(scores, "matrix")
  result_class <- class(scores)
  result <- dplyr::left_join(scores, td_summary, by = "cluster")

  td_confident <- !is.na(result$td_purity) &
    result$td_purity >= td_min_purity & result$td_mean_score >= td_min_score
  is_annotated <- result$decision == "annotate"
  conflict <- is_annotated & td_confident & (result$best_cell_type != result$td_call)

  result$best_cell_type <- ifelse(!is_annotated & td_confident, result$td_call, result$best_cell_type)
  result$decision <- ifelse(
    conflict, "review_conflict",
    ifelse(is_annotated, "annotate",
    ifelse(td_confident, "provisional", result$decision))
  )

  attr(result, "matrix") <- matrix_attr
  class(result) <- result_class
  result
}
