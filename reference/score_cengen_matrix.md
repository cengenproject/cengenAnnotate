# Score clusters against every candidate CeNGEN cell type

Computes a "dot-plot-mirrored coherence score" for each (cluster, cell
type) pair: how consistently the cluster's marker gene panel is both
highly expressed (coverage, mirroring dot *size*) and specific
(mirroring dot *color*) in that cell type, relative to CeNGEN's other
cell types. This is the computational equivalent of visually judging
whether a CeNGEN dot plot shows one cell type standing out.

## Usage

``` r
score_cengen_matrix(
  markers,
  reference,
  top_n = 20,
  min_pct_expr = 10,
  combine = c("geometric", "arithmetic"),
  cov_weight = 0.5,
  spec_weight = 0.5,
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
  `"geometric"` (default) or `"arithmetic"`.

- cov_weight, spec_weight:

  Weights on the coverage and specificity signals (exponents for
  `"geometric"`, linear weights for `"arithmetic"`).

## Value

A long tibble with columns `cluster`, `cell_type`, `score`,
`n_genes_used`, `n_genes_requested`, `n_genes_missing_from_reference`,
`n_genes_uninformative`.

## Details

For each gene `g` and cell type `t`:

- `cov(g,t) = 0` if `pct_expr(g,t) < min_pct_expr`, else
  `pct_expr(g,t) / 100`.

- `z(g,t) = (avg_expr(g,t) - mu_g) / sigma_g`, where `mu_g`/`sigma_g`
  are gene `g`'s mean/sd of `avg_expr` across all cell types.
  `spec(g,t) = plogis(z(g,t))`.

- `s(g,t)` combines the two: for `combine = "geometric"` (default),
  `cov(g,t)^cov_weight * spec(g,t)^spec_weight` (requires both signals
  to be non-trivial simultaneously); for `"arithmetic"`,
  `cov_weight * cov(g,t) + spec_weight * spec(g,t)`.

The cluster-level score for cell type `t` is the mean of `s(g,t)` over
the cluster's usable marker genes. Genes absent from the reference, or
with near-zero variance across cell types (uninformative), are excluded
from scoring and tallied in the diagnostic columns.
