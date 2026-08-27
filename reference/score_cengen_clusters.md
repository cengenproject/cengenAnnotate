# Score and classify clusters against CeNGEN cell types

Convenience wrapper that runs
[`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md)
followed by
[`classify_cengen_matches()`](https://cengenproject.github.io/cengenAnnotate/reference/classify_cengen_matches.md),
and attaches the full cluster x cell-type score matrix to the result
(retrieve it with
[`cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/cengen_matrix.md))
for drilling into all candidate cell types beyond the top 2.

## Usage

``` r
score_cengen_clusters(
  markers,
  reference,
  top_n = 20,
  min_pct_expr = 10,
  combine = c("arithmetic", "geometric"),
  cov_weight = 0.5,
  spec_weight = 0.5,
  cutoff = 0.6,
  min_gap = 0.1,
  min_genes_used = 5,
  cluster_col = "cluster",
  gene_col = "gene"
)
```

## Arguments

- markers:

  Either raw marker output (see
  [`prepare_marker_panel()`](https://cengenproject.github.io/cengenAnnotate/reference/prepare_marker_panel.md))
  or an already-prepared panel (a data frame with `cluster`, `gene`,
  `rank`, `n_used` columns).

- reference:

  A `cengen_reference` object (see
  [`new_cengen_reference()`](https://cengenproject.github.io/cengenAnnotate/reference/new_cengen_reference.md)
  /
  [`load_cengen_reference()`](https://cengenproject.github.io/cengenAnnotate/reference/load_cengen_reference.md)).

- top_n, cluster_col, gene_col:

  Passed to
  [`prepare_marker_panel()`](https://cengenproject.github.io/cengenAnnotate/reference/prepare_marker_panel.md)
  when `markers` is not already a prepared panel.

- min_pct_expr:

  Minimum percent-expressing (0-100) for a gene to count as "covered" in
  a cell type; below this, coverage is 0.

- combine:

  How to combine the coverage and specificity signals per gene:
  `"arithmetic"` (default) or `"geometric"`.

- cov_weight, spec_weight:

  Weights on the coverage and specificity signals (exponents for
  `"geometric"`, linear weights for `"arithmetic"`).

- cutoff, min_gap, min_genes_used:

  Passed to
  [`classify_cengen_matches()`](https://cengenproject.github.io/cengenAnnotate/reference/classify_cengen_matches.md).

## Value

A `cengen_scores` tibble, one row per cluster (see
[`classify_cengen_matches()`](https://cengenproject.github.io/cengenAnnotate/reference/classify_cengen_matches.md)),
with the full score matrix attached as an attribute.
