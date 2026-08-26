# Plot a CeNGEN-style dot plot for one cluster's marker panel

Reproduces the view the CeNGEN website would show for a cluster's marker
gene panel, restricted to a shortlist of candidate cell types, using the
same signals the coherence score is built from: dot size is percent of
cells expressing (coverage), dot color is the per-gene z-score
(specificity). Useful for visually inspecting `"review_ambiguous"`
clusters, or for sanity-checking an `"annotate"` call.

## Usage

``` r
plot_cengen_dotplot(
  reference,
  markers,
  cluster,
  scores = NULL,
  cell_types = NULL,
  top_n_cell_types = 8,
  top_n = 20,
  cluster_col = "cluster",
  gene_col = "gene"
)
```

## Arguments

- reference:

  A `cengen_reference` object.

- markers:

  Raw or prepared marker panel (see
  [`score_cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_matrix.md)).

- cluster:

  The cluster to plot (matched against the `cluster` column/prepared
  panel).

- scores:

  Optional `cengen_scores` object (from
  [`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md));
  if supplied and `cell_types` is not, the top `top_n_cell_types`
  candidates for `cluster` are plotted.

- cell_types:

  Optional explicit vector of cell types to plot. Required if `scores`
  is not supplied.

- top_n_cell_types:

  Number of top-scoring cell types to plot when deriving the shortlist
  from `scores` (default 8).

- top_n, cluster_col, gene_col:

  Passed to
  [`prepare_marker_panel()`](https://cengenproject.github.io/cengenAnnotate/reference/prepare_marker_panel.md)
  when `markers` is not already a prepared panel.

## Value

A `ggplot` object.
