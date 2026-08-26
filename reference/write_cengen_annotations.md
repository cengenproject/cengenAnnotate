# Write CeNGEN annotations onto a Seurat object

Opt-in step that writes
[`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md)
results onto a Seurat object as per-cell metadata columns. Never called
implicitly by the scoring functions, and by default only writes clusters
whose `decision` is `"annotate"` — clusters sent to review stay
unannotated (`NA`) unless `include` is widened.

## Usage

``` r
write_cengen_annotations(
  object,
  scores,
  cluster_col = "seurat_clusters",
  include = "annotate",
  metadata_col = "cengen_annotation",
  score_col = "cengen_annotation_score",
  overwrite = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- scores:

  A `cengen_scores` object, as returned by
  [`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md).

- cluster_col:

  Name of the cluster identity column in `object@meta.data` (default
  `"seurat_clusters"`).

- include:

  Which `decision` value(s) to write onto the object (default
  `"annotate"` only).

- metadata_col, score_col:

  Names of the new metadata columns to add, holding the matched neuron
  type and its score respectively.

- overwrite:

  If `FALSE` (default), errors if `metadata_col` already exists in
  `object@meta.data`.

## Value

The modified `Seurat` object (not mutated in place, matching Seurat's
own convention).
