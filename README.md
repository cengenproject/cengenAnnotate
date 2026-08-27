
# cengenAnnotate

<!-- badges: start -->

[![R-CMD-check](https://github.com/cengenproject/cengenAnnotate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cengenproject/cengenAnnotate/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`cengenAnnotate` automates annotation of Seurat clusters in C. elegans
single-cell data using [CeNGEN](https://www.cengen.org) reference
expression data. It replaces the manual workflow of copying marker genes
into the CeNGEN website, reading the resulting dot plot by eye, and
judging whether a cluster’s marker panel points to a single coherent
cell type.

`cengenAnnotate` computes a “dot-plot-mirrored coherence score” for
every (cluster, candidate cell type) pair directly from CeNGEN’s
unthresholded expression data (the same average-expression and
percent-expressing values a dot plot encodes as color and size), then
classifies each cluster as confidently annotatable or in need of manual
review, based on a tunable score cutoff and the gap to the next-best
candidate.

We also provide a vignette where `cengenAnnotate` is used in combination
with Seurat’s `TransferData()` to maximize annotation accuracy and depth
— see `vignette("full-annotation-protocol")`.

## Installation

``` r
# install.packages("pak")
pak::pak("cengenproject/cengenAnnotate")
```

## Workflow

``` r
library(Seurat)
library(cengenAnnotate)

# 1. Marker genes, as usual
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

# 2. Pick a CeNGEN reference dataset (stage/sex variant)
list_cengen_datasets()
reference <- load_cengen_reference("adult_herm")

# 3. Score every cluster against every candidate cell type
result <- score_cengen_clusters(markers, reference)
result

# 4. Write confident annotations back onto the Seurat object; ambiguous /
#    marker-poor clusters are left for manual review
seurat_obj <- write_cengen_annotations(seurat_obj, result)

# 5. Visually sanity-check any cluster, especially ones sent to review
plot_cengen_dotplot(reference, markers, cluster = "3", scores = result)
```

See `vignette("cengenAnnotate")` for a complete walkthrough.

## Reference data

Reference datasets (developmental stage x sex variants of CeNGEN’s
unthresholded expression data) are published to the
[`cengen-reference-data`](https://github.com/cengenproject/cengen-reference-data)
repo as a [`pins`](https://pins.rstudio.com) board, so new datasets
become available to `list_cengen_datasets()`/`load_cengen_reference()`
without any changes to this package. See that repo’s README for the
dataset schema and how to publish a new dataset.

## Defaults are provisional

The default score `cutoff` (0.6) and `min_gap` (0.10) used by
`score_cengen_clusters()` are reasonable starting points, not
empirically calibrated thresholds. Validate them against clusters with
independently known identity before relying on this package for
unattended annotation.
