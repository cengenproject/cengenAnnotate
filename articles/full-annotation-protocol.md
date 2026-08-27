# Protocol: annotating a new dataset end to end

The [“Introduction to
cengenAnnotate”](https://cengenproject.github.io/cengenAnnotate/articles/cengenAnnotate.md)
vignette walks through the package’s functions in isolation on synthetic
data. This vignette is different: it’s the **protocol**, in order, for
taking a new C. elegans single-cell dataset from raw counts to a Seurat
object with reviewed cell-type annotations. It’s exactly the sequence of
steps used (twice, independently) to annotate an L1 and an adult
hermaphrodite whole-animal dataset — including the steps that aren’t
`cengenAnnotate` function calls at all, because getting those right is
what makes the scoring step work well.

Every code chunk here is `eval = FALSE` — it’s meant to be adapted to
your own object names and copied into your own script, not run as-is.

## Overview

1.  Standard Seurat QC and preprocessing
2.  Cluster and visualize
3.  *(Optional, recommended for whole-animal data)* Subset to the
    lineage of interest
4.  Find markers and score with `cengenAnnotate`
5.  Resolve ambiguous clusters by targeted reclustering
6.  Cross-validate with
    [`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html)
7.  Combine into a 3-tier call and review what’s left
8.  Write annotations back onto the Seurat object

Steps 4–8 are what `cengenAnnotate` automates or supports directly.
Steps 1–3 are standard Seurat work, included here because the quality of
everything downstream depends on getting them right for *your* data.

## 1. QC and preprocessing

Standard Seurat QC: filter on `nFeature_RNA`, `nCount_RNA`, and
`percent.mt` (or the mitochondrial-equivalent for your prep), using
thresholds appropriate to your sequencing depth and tissue.

``` r

obj <- Seurat::CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
```

**Doublet detection is optional QC hygiene, not a validated fix for
ambiguous clusters.** Some of the ambiguous “blob” clusters described in
step 5 do have expression profiles consistent with doublets (elevated
`nCount`/`nFeature`, lower `percent.mt`, concentrated in a handful of
specific clusters rather than spread evenly) — but when this was tested
directly, removing those cells and re-scoring left 14 of 15 affected
clusters’ `cengenAnnotate` calls completely unchanged, and made the 15th
slightly worse (most likely just a smaller-sample-size effect, not a
real one). So: worth running for its own sake, but don’t expect it to
resolve ambiguous calls on its own.

``` r

sce <- scDblFinder::scDblFinder(Seurat::GetAssayData(obj, layer = "counts"))
obj$is_doublet <- sce$scDblFinder.class == "doublet"
obj <- subset(obj, subset = !is_doublet)
```

## 2. Normalize, scale, cluster, visualize

[`Seurat::SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)
is the modern replacement for the
[`NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
/
[`FindVariableFeatures()`](https://satijalab.org/seurat/reference/FindVariableFeatures.html)
/ [`ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
sequence — it regresses out sequencing depth as part of its own
regularized model, and takes `percent.mt` as an explicit covariate to
regress out too (depth alone doesn’t remove a mitochondrial-fraction
effect):

``` r

obj <- Seurat::SCTransform(obj, vars.to.regress = "percent.mt")
obj <- Seurat::RunPCA(obj)
obj <- Seurat::FindNeighbors(obj, dims = 1:30)

# algorithm = 1 (Louvain, default) or algorithm = 4 (Leiden, requires
# reticulate + the Python leidenalg/python-igraph packages). Leiden
# guarantees well-connected communities where Louvain can occasionally
# produce disconnected ones; in practice, on both reference datasets the
# two produced very similar downstream cengenAnnotate results, so this
# is not a decision to agonize over.
obj <- Seurat::FindClusters(obj, algorithm = 1)
obj <- Seurat::RunUMAP(obj, dims = 1:30)
Seurat::DimPlot(obj, label = TRUE)
```

**A caveat on switching to
[`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html):**
every threshold elsewhere in this protocol — the 65%-coverage
neuron-panel cutoff in step 3, the clustering resolution target in step
4, and `cengenAnnotate`’s own `cutoff`/`min_gap` defaults — was tuned
against the classic `NormalizeData`/`ScaleData` pipeline, not against
[`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)
output. Switching normalization methods can shift cluster boundaries and
marker gene rankings, so re-check those thresholds (the same way
described in step 3) rather than assuming they carry over unchanged.

Look at the UMAP before going further. Clusters that are visually
diffuse, span multiple separated islands, or sit in a long “smear”
between two otherwise-distinct clusters are worth noting now — they’re
the ones most likely to need the targeted reclustering in step 5.

## 3. Subset to the lineage of interest (optional)

If you only care about one class of cells (this protocol was developed
annotating neurons out of whole-animal data), don’t try to score every
cluster in the whole dataset against a neuron-only reference. Subset
first — and do the classification at the **cluster** level, using the
clusters from step 2, not per-cell.

The pitfall: picking a single marker gene to threshold on. **No single
gene is expressed in every cell type in the class you’re after** — pick
a gene that’s strong in most neuron types and you will still cleanly
lose whichever types happen to express it weakly. A panel fixes this,
but only if it’s built to fix exactly that failure mode: for each
candidate gene, look at its *worst* covered target cell type, not its
average. Build the panel by greedy set-cover against the CeNGEN
reference — repeatedly find the target cell type the current panel
covers worst, and add whichever candidate gene rescues that specific
type the most, until coverage stops improving:

``` r

# `reference` is the same cengen_reference used for scoring in step 5.
# `target_types` is every CeNGEN cell type in the class you're
# subsetting to (e.g. every neuron type); `candidate_pool` is every gene
# worth considering (e.g. neuropeptide/synaptic genes for neurons) --
# restrict this so the greedy search doesn't just rediscover one
# ubiquitously-expressed housekeeping gene.
coverage_for_panel <- function(panel) {
  reference$data |>
    dplyr::filter(gene %in% panel, cell_type %in% target_types) |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(coverage = max(pct_expr), .groups = "drop")
}

panel <- character(0)
repeat {
  cov <- coverage_for_panel(panel)
  worst <- cov[which.min(cov$coverage), ]
  # every target type already well covered -- or no candidate improves
  # the worst point any further -- stop
  candidates <- reference$data |>
    dplyr::filter(gene %in% candidate_pool, !gene %in% panel, cell_type == worst$cell_type) |>
    dplyr::arrange(dplyr::desc(pct_expr))
  if (nrow(cov) > 0 && nrow(candidates) > 0 && candidates$pct_expr[1] > worst$coverage) {
    panel <- c(panel, candidates$gene[1])
  } else {
    break
  }
}
```

Then classify each existing cluster as in-class if *any* panel gene is
detected in a high enough fraction of that cluster’s cells (OR logic
across genes, applied per cluster — not a per-cell “expresses at least
one panel gene” filter, which is far noisier given single-cell dropout):

``` r

panel_present <- intersect(panel, rownames(obj))
pos_mat <- Seurat::FetchData(obj, vars = panel_present, layer = "counts") > 0
pct_by_cluster <- sapply(panel_present, function(g) {
  tapply(pos_mat[, g], obj$seurat_clusters, mean) * 100
})
max_pct <- apply(pct_by_cluster, 1, max)

target_clusters <- names(max_pct)[max_pct >= 65]  # see below for where 65 comes from
neurons <- subset(obj, subset = seurat_clusters %in% target_clusters)
```

Before trusting a threshold like the 65% above, check it against the
reference: compute each *known* target cell type’s best coverage under
the finished panel (the same `coverage_for_panel()` from the build step)
and each known non-target cell type’s best coverage, and confirm there’s
a sharp break between the lowest-covered true-positive type and the
highest-covered true-negative type. Pick the threshold in that gap. If
the break is sharp, a simple threshold is trustworthy; if it’s gradual,
the panel needs another gene (go back to the greedy step), not a fiddled
threshold. The 22-gene panel this protocol validated for C. elegans
neurons (built from CeNGEN’s neuropeptide/synaptic-gene reference data)
was
`sbt-1, K02F2.5, T27C4.1, Y44A6D.2, pghm-1, snt-4, ric-4, cpx-1, egl-21, cla-1, Y23H5B.1, Y39A1A.27, unc-129, nlp-3, nlp-62, nlp-58, R102.2, acbp-6, twk-49, F14B6.2, ntc-1, R74.10`,
with a 65%-coverage threshold — reusable as a starting point for other
C. elegans datasets, but still worth re-checking the gap on your own
data rather than trusted blindly.

## 4. Recluster the subset

Whatever clustering resolution and PCs were used on the whole dataset
are almost never right for a subset — a lineage that was 3 clusters in
the full dataset usually needs finer resolution once it’s the only thing
in the embedding. Redo variable features, scaling, PCA, and clustering
on the subset from scratch, typically at higher `resolution`:

``` r

neurons <- Seurat::SCTransform(neurons, vars.to.regress = "percent.mt")
neurons <- Seurat::RunPCA(neurons)
neurons <- Seurat::FindNeighbors(neurons, dims = 1:30)
neurons <- Seurat::FindClusters(neurons, resolution = 4)
neurons <- Seurat::RunUMAP(neurons, dims = 1:30)
```

## 5. Score clusters against CeNGEN

This is the part `cengenAnnotate` automates directly.

``` r

library(cengenAnnotate)

# SCT-normalized counts need PrepSCTFindMarkers() first, so fold changes
# are comparable across cells that got different depth corrections --
# skip this line if you stuck with the classic NormalizeData()/ScaleData()
# pipeline instead of SCTransform() in steps 2/4.
neurons <- Seurat::PrepSCTFindMarkers(neurons)
markers <- Seurat::FindAllMarkers(neurons, assay = "SCT", only.pos = TRUE)
reference <- load_cengen_reference("L4_herm")  # match developmental stage/sex where possible

result <- score_cengen_clusters(markers, reference)
result  # sorted by confidence, not cluster number
```

`result$decision` is one of `"annotate"`, `"review_ambiguous"`, or
`"review_insufficient_markers"`. Inspect the ones sent to review with
the full cluster x cell-type matrix and the QC dot plot before deciding
what to do with them:

``` r

cengen_matrix(result) |> dplyr::filter(cluster == "12") |> dplyr::arrange(dplyr::desc(score))
plot_cengen_dotplot(reference, markers, cluster = "12", scores = result)
```

## 6. Resolve ambiguous clusters by targeted reclustering

Some `review_ambiguous` clusters resolve cleanly if you subset just that
cluster (or a small group of adjacent, similarly-ambiguous clusters —
“blobs” that touch or overlap on the UMAP) and recluster it alone at
higher resolution, then re-score just the new sub-clusters against the
reference.

``` r

blob <- subset(neurons, idents = c("12", "19", "31"))
blob <- Seurat::FindVariableFeatures(blob)
blob <- Seurat::ScaleData(blob, features = Seurat::VariableFeatures(blob))
blob <- Seurat::RunPCA(blob)
blob <- Seurat::FindNeighbors(blob, dims = 1:20)
blob <- Seurat::FindClusters(blob, resolution = 6)

blob_markers <- Seurat::FindAllMarkers(blob, only.pos = TRUE)
blob_result <- score_cengen_clusters(blob_markers, reference)
```

Not every ambiguous cluster resolves this way, and it’s worth knowing
the two different reasons one won’t, since they call for different
responses:

- **Genuine transcriptional continua.** Some ventral cord motor neuron
  families (e.g. VA/VB/DA/DB/AS) blend into each other transcriptionally
  regardless of clustering algorithm or resolution — this is real
  biology, not a pipeline failure, and no amount of reclustering will
  produce a clean split. These are expected to stay ambiguous.
- **Heterogeneous mixtures with normal QC metrics.** A cluster that
  mixes clearly unrelated cell types, but whose cells look fine on
  standard QC, is a signal worth investigating further (see the
  doublet-detection note in step 1) rather than something to keep
  reclustering around.

## 7. Cross-validate with Seurat’s TransferData

`cengenAnnotate`’s own scoring is conservative by design: it strongly
prefers leaving a cluster unannotated over guessing wrong. That means it
under-covers relative to
[`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html),
which makes a confident-looking call for nearly every cluster but is
less precise on the calls it’s actually unsure about. Running both and
comparing is more informative than trusting either alone.

``` r

ref_seurat <- load_reference_seurat_object("L4_herm")  # your own reference Seurat object, cell_type in metadata
anchors <- Seurat::FindTransferAnchors(reference = ref_seurat, query = neurons)
predictions <- Seurat::TransferData(anchorset = anchors, refdata = ref_seurat$cell_type)
neurons <- Seurat::AddMetaData(neurons, metadata = predictions)
```

## 8. Combine into a 3-tier call

[`combine_with_transfer_data()`](https://cengenproject.github.io/cengenAnnotate/reference/combine_with_transfer_data.md)
merges the two into three tiers: `cengenAnnotate`’s own `"annotate"`
calls are left untouched; clusters it left ambiguous but where
`TransferData`’s per-cell vote is confident (majority-vote purity and
mean prediction score both above threshold) are promoted to
`"provisional"`, labeled with `TransferData`’s call; a genuine conflict
— both methods confident, different answers — is flagged
`"review_conflict"` rather than silently resolved one way or the other.

``` r

combined <- combine_with_transfer_data(result, neurons)
table(combined$decision)
```

On both datasets this protocol was validated against, `"annotate"` calls
were essentially always correct (99–100%), `"provisional"` calls were
good but noticeably less precise (81–93%), and the two together covered
82–87% of all clusters — leaving the rest genuinely unresolved. Treat
those percentages as a sense of what to expect, not a guarantee:
re-check precision on your own data if you have any independent way to
(a known subset of markers, a published reference for the same tissue,
etc.), since both numbers were measured against paper-published ground
truth that won’t exist for every new dataset.

## 9. Manual review of what’s left

Whatever remains `"review_ambiguous"`, `"review_insufficient_markers"`,
or `"review_conflict"` after steps 6–8 needs a human.
[`cengen_matrix()`](https://cengenproject.github.io/cengenAnnotate/reference/cengen_matrix.md)
and
[`plot_cengen_dotplot()`](https://cengenproject.github.io/cengenAnnotate/reference/plot_cengen_dotplot.md)
are the same tools used for triage during scoring in step 5 — use them
here on the clusters that made it this far without a confident call from
either method.

## 10. Write annotations back

[`write_cengen_annotations()`](https://cengenproject.github.io/cengenAnnotate/reference/write_cengen_annotations.md)
is opt-in per tier — nothing is written for clusters you don’t
explicitly `include`:

``` r

neurons <- write_cengen_annotations(neurons, combined, include = c("annotate", "provisional"))
table(neurons$cengen_annotation, useNA = "ifany")
```

Consider keeping `"annotate"` and `"provisional"` clusters visibly
distinguished downstream (e.g. a separate metadata column, or carrying
`combined$decision` through) rather than merging them into one
undifferentiated label — the accuracy gap between the two tiers means
they carry different amounts of trust.

## A note on cutoffs

The default `cutoff`, `min_gap`, and `combine = "arithmetic"` in
[`score_cengen_clusters()`](https://cengenproject.github.io/cengenAnnotate/reference/score_cengen_clusters.md)
reflect what worked on the two datasets this protocol was validated on,
not a universal calibration. Similarly, the
`td_min_purity`/`td_min_score` defaults in
[`combine_with_transfer_data()`](https://cengenproject.github.io/cengenAnnotate/reference/combine_with_transfer_data.md)
were chosen empirically for those same runs. If you have any independent
way to check precision on your own data, do — and adjust the cutoffs
before trusting the pipeline unattended on a new stage, sex, or tissue
context that hasn’t been checked this way before.
