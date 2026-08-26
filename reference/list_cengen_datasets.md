# List available CeNGEN reference datasets

Lists the CeNGEN reference datasets currently published on the pins
board, along with their metadata (developmental stage, sex, when they
were prepared, etc.), without downloading the (potentially large)
underlying expression tables.

## Usage

``` r
list_cengen_datasets(board = cengen_board(), pattern = NULL)
```

## Arguments

- board:

  A pins board, defaulting to
  [`cengen_board()`](https://cengenproject.github.io/cengenAnnotate/reference/cengen_board.md).

- pattern:

  Optional regular expression to filter dataset names.

## Value

A tibble with one row per dataset: `name`, `stage`, `sex`,
`dataset_label`, `source_version`, `date_prepared`, `n_cell_types`,
`n_genes`, `version`, `created`.
