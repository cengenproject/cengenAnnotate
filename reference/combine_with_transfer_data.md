# Combine cengenAnnotate results with Seurat TransferData predictions

Augments a `cengen_scores` object (from
[`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md))
with per-cluster predictions from
[`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html),
producing a 3-tier classification validated across two independent
datasets: clusters `cengenAnnotate` already confidently called stay
`"annotate"`; clusters it didn't, but where TransferData's own per-cell
votes are confident, get promoted to `"provisional"` (labeled with
TransferData's call); everything else is left as-is for manual review.

## Usage

``` r
combine_with_transfer_data(
  scores,
  object,
  cluster_col = "seurat_clusters",
  predicted_col = "predicted.id",
  score_col = "prediction.score.max",
  td_min_purity = 0.6,
  td_min_score = 0.6
)
```

## Arguments

- scores:

  A `cengen_scores` object, as returned by
  [`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md).

- object:

  A Seurat object with TransferData results already added (i.e. has
  `predicted_col` and `score_col` metadata columns).

- cluster_col:

  Name of the cluster identity column in `object@meta.data` (default
  `"seurat_clusters"`).

- predicted_col, score_col:

  Names of the per-cell TransferData columns (Seurat's own defaults:
  `"predicted.id"`, `"prediction.score.max"`).

- td_min_purity:

  Minimum fraction of a cluster's cells that must agree on
  TransferData's majority vote for that cluster to count as
  TransferData-confident.

- td_min_score:

  Minimum mean per-cell `prediction.score.max` within a cluster for it
  to count as TransferData-confident.

## Value

A `cengen_scores` tibble (same shape as `scores`, with the score matrix
attribute preserved so
[`cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/cengen_matrix.md)
still works), with `decision` extended to include `"provisional"` and
`"review_conflict"`, `best_cell_type` updated to TransferData's call for
`"provisional"` rows, and three new diagnostic columns: `td_call`,
`td_purity`, `td_mean_score`.

## Details

On both test datasets, `cengenAnnotate`'s `"annotate"` calls and
TransferData's confident per-cell-vote majority never disagreed. A
genuine conflict (both methods confident, different answers) is
therefore treated as a rare edge case rather than silently resolved:
it's flagged as `"review_conflict"`, keeping `cengenAnnotate`'s own
(independently validated) call as `best_cell_type` while recording
TransferData's call in `td_call` for inspection.

`object` must already have
[`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html)
results attached, e.g.:

    anchors <- Seurat::FindTransferAnchors(reference = ref, query = object)
    predictions <- Seurat::TransferData(anchorset = anchors, refdata = ref$cell_type)
    object <- Seurat::AddMetaData(object, metadata = predictions)
