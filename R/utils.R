`%||%` <- function(x, y) if (is.null(x)) y else x

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
