#' List available CeNGEN reference datasets
#'
#' Lists the CeNGEN reference datasets currently published on the pins
#' board, along with their metadata (developmental stage, sex, when they
#' were prepared, etc.), without downloading the (potentially large)
#' underlying expression tables.
#'
#' @param board A pins board, defaulting to [cengen_board()].
#' @param pattern Optional regular expression to filter dataset names.
#'
#' @return A tibble with one row per dataset: `name`, `stage`, `sex`,
#'   `dataset_label`, `source_version`, `date_prepared`, `n_cell_types`,
#'   `n_genes`, `version`, `created`.
#' @export
list_cengen_datasets <- function(board = cengen_board(), pattern = NULL) {
  names <- pins::pin_list(board)
  if (!is.null(pattern)) {
    names <- names[grepl(pattern, names)]
  }

  if (length(names) == 0) {
    return(tibble::tibble(
      name = character(), stage = character(), sex = character(),
      dataset_label = character(), source_version = character(),
      date_prepared = character(), n_cell_types = integer(),
      n_genes = integer(), version = character(), created = as.POSIXct(character())
    ))
  }

  rows <- lapply(names, function(nm) {
    meta <- pins::pin_meta(board, nm)
    user <- meta$user %||% list()
    tibble::tibble(
      name = nm,
      stage = as.character(user$stage %||% NA),
      sex = as.character(user$sex %||% NA),
      dataset_label = as.character(user$dataset_label %||% NA),
      source_version = as.character(user$source_version %||% NA),
      date_prepared = as.character(user$date_prepared %||% NA),
      n_cell_types = as.integer(user$n_cell_types %||% NA),
      n_genes = as.integer(user$n_genes %||% NA),
      version = as.character(meta$version$version %||% NA),
      created = meta$created %||% as.POSIXct(NA)
    )
  })

  dplyr::bind_rows(rows)
}
