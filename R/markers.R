#' Prepare a per-cluster marker gene panel
#'
#' Takes raw [Seurat::FindAllMarkers()] output (or any per-cluster marker
#' `data.frame` with the same shape) and reduces it to the top marker genes
#' per cluster used for CeNGEN scoring. Fold-change and adjusted p-value are
#' used only to *select* which genes make the panel; once selected, every
#' gene in the panel is treated equally by [score_cengen_matrix()].
#'
#' @param markers A data frame of per-cluster markers, e.g. the output of
#'   `Seurat::FindAllMarkers()`.
#' @param top_n Number of top marker genes to keep per cluster (default 20).
#' @param cluster_col,gene_col Column names identifying the cluster and gene.
#' @param p_val_adj_col,avg_log2FC_col Column names used for filtering/ranking.
#'   Set to `NULL` to skip the corresponding filter/ranking step.
#' @param max_p_val_adj Maximum adjusted p-value to keep. Ignored if
#'   `p_val_adj_col` is `NULL`.
#' @param only_pos If `TRUE` (default), keep only genes with positive
#'   `avg_log2FC` (upregulated markers). Ignored if `avg_log2FC_col` is `NULL`.
#'
#' @return A tibble with columns `cluster`, `gene`, `rank`, `n_requested`,
#'   `n_used` (the last two constant within a cluster; `n_used` is the
#'   number of marker genes actually kept for that cluster, which may be
#'   less than `top_n` if fewer significant markers were available).
#' @export
prepare_marker_panel <- function(
  markers,
  top_n = 20,
  cluster_col = "cluster",
  gene_col = "gene",
  p_val_adj_col = "p_val_adj",
  avg_log2FC_col = "avg_log2FC",
  max_p_val_adj = 0.05,
  only_pos = TRUE
) {
  required_cols <- c(cluster_col, gene_col)
  if (!is.null(p_val_adj_col) && !is.null(max_p_val_adj)) {
    required_cols <- c(required_cols, p_val_adj_col)
  }
  if (!is.null(avg_log2FC_col) && only_pos) {
    required_cols <- c(required_cols, avg_log2FC_col)
  }
  assert_has_columns(markers, required_cols, "markers", "cengen_invalid_markers_error")

  markers <- tibble::as_tibble(markers)
  markers$.cluster <- as.character(markers[[cluster_col]])
  markers$.gene <- normalize_gene_symbol(as.character(markers[[gene_col]]))

  if (!is.null(avg_log2FC_col) && only_pos) {
    markers <- markers[markers[[avg_log2FC_col]] > 0, , drop = FALSE]
  }
  if (!is.null(p_val_adj_col) && !is.null(max_p_val_adj)) {
    markers <- markers[markers[[p_val_adj_col]] <= max_p_val_adj, , drop = FALSE]
  }

  sort_cols <- list()
  if (!is.null(avg_log2FC_col)) sort_cols$fc <- -markers[[avg_log2FC_col]]
  if (!is.null(p_val_adj_col)) sort_cols$p <- markers[[p_val_adj_col]]
  if (length(sort_cols) > 0) {
    markers <- markers[do.call(order, unname(sort_cols)), , drop = FALSE]
  }

  split_by_cluster <- split(markers, markers$.cluster)
  panel <- lapply(names(split_by_cluster), function(cl) {
    rows <- split_by_cluster[[cl]]
    rows <- rows[!duplicated(rows$.gene), , drop = FALSE]
    n_used <- min(nrow(rows), top_n)
    rows <- utils::head(rows, top_n)
    tibble::tibble(
      cluster = cl,
      gene = rows$.gene,
      rank = seq_len(nrow(rows)),
      n_requested = top_n,
      n_used = n_used
    )
  })

  dplyr::bind_rows(panel)
}

is_prepared_panel <- function(markers) {
  all(c("cluster", "gene", "rank", "n_used") %in% names(markers))
}
