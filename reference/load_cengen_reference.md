# Load a CeNGEN reference dataset from the pins board

Fetches a named dataset (e.g. `"adult_herm"`) from a CeNGEN `pins` board
and wraps it as a
[`new_cengen_reference()`](https://cengenproject.github.io/cengenAnnotate/reference/new_cengen_reference.md)
object. Remote boards are cached locally by `pins`, so repeated calls
for the same version do not re-download.

## Usage

``` r
load_cengen_reference(name, board = cengen_board(), version = NULL)
```

## Arguments

- name:

  Dataset name, as returned by
  [`list_cengen_datasets()`](https://cengenproject.github.io/cengenAnnotate/reference/list_cengen_datasets.md).

- board:

  A pins board, defaulting to
  [`cengen_board()`](https://cengenproject.github.io/cengenAnnotate/reference/cengen_board.md).

- version:

  Optional specific pin version; defaults to the latest.

## Value

A `cengen_reference` object.
