# Classify clusters as annotatable or in need of review

Applies the cutoff/gap decision rule to a cluster x cell-type score
matrix (from
[`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md)),
reducing each cluster to its best and second-best matching cell type and
a decision.

## Usage

``` r
classify_cengen_matches(
  score_matrix,
  cutoff = 0.6,
  min_gap = 0.1,
  min_genes_used = 5
)
```

## Arguments

- score_matrix:

  A long score matrix, as returned by
  [`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md).

- cutoff:

  Minimum `best_score` required to auto-annotate.

- min_gap:

  Minimum gap between the best and second-best score required to
  auto-annotate (guards against ambiguous ties).

- min_genes_used:

  Minimum number of usable marker genes required; below this, the
  cluster is always sent to review.

## Value

A tibble with one row per cluster: `cluster`, `best_cell_type`,
`best_score`, `second_cell_type`, `second_score`, `gap`, `n_genes_used`,
`n_genes_requested`, `n_genes_missing_from_reference`,
`n_genes_uninformative`, `decision`. Sorted by `best_score` descending
(most confident cluster first), not by input cluster order.

## Details

Decision priority, evaluated in order:

1.  If `n_genes_used < min_genes_used`: `"review_insufficient_markers"`,
    regardless of score (too few usable marker genes to trust any
    match).

2.  Else if `best_score >= cutoff` AND `gap >= min_gap`: `"annotate"`
    (one cell type stands out clearly).

3.  Otherwise: `"review_ambiguous"` (either no candidate scores high
    enough, or more than one candidate is plausible).
