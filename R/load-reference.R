#' Construct a CeNGEN reference object
#'
#' Wraps a tidy long table of unthresholded CeNGEN expression data (as
#' published on the `cengen-reference-data` pins board) into a
#' `cengen_reference` object, precomputing the per-gene statistics that
#' [score_cengen_matrix()] needs.
#'
#' @param data A data frame with columns `gene`, `cell_type`, `avg_expr`
#'   (unthresholded average expression level) and `pct_expr` (percent of
#'   cells of that cell type expressing the gene, on a 0-100 scale).
#' @param meta A named list of metadata describing the dataset (e.g.
#'   `stage`, `sex`, `source_version`). Optional.
#'
#' @return A `cengen_reference` object.
#' @export
new_cengen_reference <- function(data, meta = list()) {
  assert_has_columns(
    data, c("gene", "cell_type", "avg_expr", "pct_expr"),
    "data", "cengen_invalid_reference_error"
  )

  data <- tibble::as_tibble(data)
  data$gene <- normalize_gene_symbol(as.character(data$gene))
  data$cell_type <- as.character(data$cell_type)
  data$avg_expr <- as.numeric(data$avg_expr)
  data$pct_expr <- as.numeric(data$pct_expr)

  dup_key <- paste(data$gene, data$cell_type, sep = "\r")
  if (anyDuplicated(dup_key) > 0) {
    cengen_abort(
      "`data` has duplicate (gene, cell_type) rows; each pair must appear at most once.",
      class = "cengen_invalid_reference_error"
    )
  }

  gene_stats <- dplyr::summarise(
    dplyr::group_by(data, .data$gene),
    mu = mean(.data$avg_expr, na.rm = TRUE),
    sigma = stats::sd(.data$avg_expr, na.rm = TRUE),
    .groups = "drop"
  )

  structure(
    list(data = data, meta = meta, gene_stats = gene_stats),
    class = "cengen_reference"
  )
}

#' @export
print.cengen_reference <- function(x, ...) {
  meta <- x$meta
  cat("<cengen_reference>\n")
  cat(sprintf("  stage:   %s\n", meta$stage %||% "?"))
  cat(sprintf("  sex:     %s\n", meta$sex %||% "?"))
  cat(sprintf("  label:   %s\n", meta$dataset_label %||% "?"))
  cat(sprintf("  prepared: %s\n", meta$date_prepared %||% "?"))
  cat(sprintf(
    "  genes: %d, cell types: %d\n",
    dplyr::n_distinct(x$data$gene), dplyr::n_distinct(x$data$cell_type)
  ))
  invisible(x)
}

#' Load a CeNGEN reference dataset from the pins board
#'
#' Fetches a named dataset (e.g. `"adult_herm"`) from a CeNGEN `pins` board
#' and wraps it as a [new_cengen_reference()] object. Remote boards are
#' cached locally by `pins`, so repeated calls for the same version do not
#' re-download.
#'
#' @param name Dataset name, as returned by [list_cengen_datasets()].
#' @param board A pins board, defaulting to [cengen_board()].
#' @param version Optional specific pin version; defaults to the latest.
#'
#' @return A `cengen_reference` object.
#' @export
load_cengen_reference <- function(name, board = cengen_board(), version = NULL) {
  data <- pins::pin_read(board, name, version = version)
  meta <- tryCatch(
    pins::pin_meta(board, name, version = version)$user,
    error = function(e) list()
  )
  new_cengen_reference(data, meta = meta %||% list())
}
