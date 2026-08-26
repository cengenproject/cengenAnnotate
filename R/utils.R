`%||%` <- function(x, y) if (is.null(x)) y else x

# Seurat replaces underscores in feature names with dashes (silently, with a
# warning) when building a Seurat object, so a gene symbol that contains an
# underscore in its "canonical" form (e.g. from a reference table) would
# otherwise fail to match the same gene as it appears in FindAllMarkers()
# output. Normalizing both sides to dashes here avoids that silent mismatch.
normalize_gene_symbol <- function(x) gsub("_", "-", x, fixed = TRUE)

cengen_abort <- function(message, class, ...) {
  rlang::abort(message, class = c(class, "cengen_error"), ...)
}

assert_has_columns <- function(data, columns, arg_name, error_class) {
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols) > 0) {
    cengen_abort(
      sprintf(
        "`%s` is missing required column(s): %s",
        arg_name, paste(missing_cols, collapse = ", ")
      ),
      class = error_class
    )
  }
  invisible(TRUE)
}
