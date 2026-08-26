#' Write CeNGEN annotations onto a Seurat object
#'
#' Opt-in step that writes [score_cengen_clusters()] results onto a Seurat
#' object as per-cell metadata columns. Never called implicitly by the
#' scoring functions, and by default only writes clusters whose `decision`
#' is `"annotate"` — clusters sent to review stay unannotated (`NA`) unless
#' `include` is widened.
#'
#' @param object A `Seurat` object.
#' @param scores A `cengen_scores` object, as returned by
#'   [score_cengen_clusters()].
#' @param cluster_col Name of the cluster identity column in
#'   `object@meta.data` (default `"seurat_clusters"`).
#' @param include Which `decision` value(s) to write onto the object
#'   (default `"annotate"` only).
#' @param metadata_col,score_col Names of the new metadata columns to add,
#'   holding the matched neuron type and its score respectively.
#' @param overwrite If `FALSE` (default), errors if `metadata_col` already
#'   exists in `object@meta.data`.
#'
#' @return The modified `Seurat` object (not mutated in place, matching
#'   Seurat's own convention).
#' @export
write_cengen_annotations <- function(
  object,
  scores,
  cluster_col = "seurat_clusters",
  include = "annotate",
  metadata_col = "cengen_annotation",
  score_col = "cengen_annotation_score",
  overwrite = FALSE
) {
  if (!inherits(object, "Seurat")) {
    cengen_abort("`object` must be a Seurat object.", class = "cengen_invalid_object_error")
  }
  assert_has_columns(object@meta.data, cluster_col, "object@meta.data", "cengen_invalid_object_error")
  assert_has_columns(
    scores, c("cluster", "best_neuron_type", "best_score", "decision"),
    "scores", "cengen_invalid_scores_error"
  )

  if (metadata_col %in% colnames(object@meta.data) && !overwrite) {
    cengen_abort(
      sprintf("Metadata column '%s' already exists; pass overwrite = TRUE to replace it.", metadata_col),
      class = "cengen_metadata_exists_error"
    )
  }

  cluster_ids <- as.character(object@meta.data[[cluster_col]])
  score_clusters <- as.character(scores$cluster)

  unmatched_in_object <- setdiff(unique(cluster_ids), score_clusters)
  if (length(unmatched_in_object) > 0) {
    warning(sprintf(
      "%d cluster(s) present in `object` have no matching row in `scores`: %s",
      length(unmatched_in_object), paste(unmatched_in_object, collapse = ", ")
    ))
  }
  unmatched_in_scores <- setdiff(score_clusters, unique(cluster_ids))
  if (length(unmatched_in_scores) > 0) {
    warning(sprintf(
      "%d cluster(s) in `scores` have no matching cluster in `object`: %s",
      length(unmatched_in_scores), paste(unmatched_in_scores, collapse = ", ")
    ))
  }

  eligible <- scores[scores$decision %in% include, , drop = FALSE]
  lookup_annotation <- stats::setNames(eligible$best_neuron_type, eligible$cluster)
  lookup_score <- stats::setNames(eligible$best_score, eligible$cluster)

  annotation <- unname(lookup_annotation[cluster_ids])
  score_vals <- unname(lookup_score[cluster_ids])

  object <- SeuratObject::AddMetaData(object, metadata = annotation, col.name = metadata_col)
  object <- SeuratObject::AddMetaData(object, metadata = score_vals, col.name = score_col)

  object
}
