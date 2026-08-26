# Construct a CeNGEN reference object

Wraps a tidy long table of unthresholded CeNGEN expression data (as
published on the `cengen-reference-data` pins board) into a
`cengen_reference` object, precomputing the per-gene statistics that
[`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md)
needs.

## Usage

``` r
new_cengen_reference(data, meta = list())
```

## Arguments

- data:

  A data frame with columns `gene`, `cell_type`, `avg_expr`
  (unthresholded average expression level) and `pct_expr` (percent of
  cells of that cell type expressing the gene, on a 0-100 scale).

- meta:

  A named list of metadata describing the dataset (e.g. `stage`, `sex`,
  `source_version`). Optional.

## Value

A `cengen_reference` object.
