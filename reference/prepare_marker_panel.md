# Prepare a per-cluster marker gene panel

Takes raw
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
output (or any per-cluster marker `data.frame` with the same shape) and
reduces it to the top marker genes per cluster used for CeNGEN scoring.
Fold-change and adjusted p-value are used only to *select* which genes
make the panel; once selected, every gene in the panel is treated
equally by
[`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md).

## Usage

``` r
prepare_marker_panel(
  markers,
  top_n = 20,
  cluster_col = "cluster",
  gene_col = "gene",
  p_val_adj_col = "p_val_adj",
  avg_log2FC_col = "avg_log2FC",
  max_p_val_adj = 0.05,
  only_pos = TRUE
)
```

## Arguments

- markers:

  A data frame of per-cluster markers, e.g. the output of
  [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html).

- top_n:

  Number of top marker genes to keep per cluster (default 20).

- cluster_col, gene_col:

  Column names identifying the cluster and gene.

- p_val_adj_col, avg_log2FC_col:

  Column names used for filtering/ranking. Set to `NULL` to skip the
  corresponding filter/ranking step.

- max_p_val_adj:

  Maximum adjusted p-value to keep. Ignored if `p_val_adj_col` is
  `NULL`.

- only_pos:

  If `TRUE` (default), keep only genes with positive `avg_log2FC`
  (upregulated markers). Ignored if `avg_log2FC_col` is `NULL`.

## Value

A tibble with columns `cluster`, `gene`, `rank`, `n_requested`, `n_used`
(the last two constant within a cluster; `n_used` is the number of
marker genes actually kept for that cluster, which may be less than
`top_n` if fewer significant markers were available).
