# Get the full cluster x cell-type score matrix

Retrieves the score matrix attached to a `cengen_scores` object by
[`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md).

## Usage

``` r
cengen_matrix(x)
```

## Arguments

- x:

  A `cengen_scores` object.

## Value

The attached long score matrix (see
[`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md)),
or `NULL` if none is attached.
