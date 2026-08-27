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

**Doublet detection is not part of this protocol.** It’s a reasonable
thing to add — some of the ambiguous “blob” clusters described in step 5
have expression profiles consistent with doublets — but it wasn’t
validated as part of this pipeline. If you skip it (as both reference
runs of this protocol did), be aware that heterogeneous clusters with
otherwise-normal QC metrics may be doublets rather than genuine
ambiguity, and a doublet-detection step (e.g. `scDblFinder`,
`DoubletFinder`) run before clustering would be the natural place to add
it if it becomes a problem for your data.

## 2. Normalize, scale, cluster, visualize

``` r

obj <- Seurat::NormalizeData(obj)
obj <- Seurat::FindVariableFeatures(obj)
# Scale only the HVGs, not all genes -- scaling every gene on a
# whole-animal object is a common cause of memory blowups and buys
# nothing downstream, since PCA only uses the HVGs anyway.
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj))
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

Look at the UMAP before going further. Clusters that are visually
diffuse, span multiple separated islands, or sit in a long “smear”
between two otherwise-distinct clusters are worth noting now — they’re
the ones most likely to need the targeted reclustering in step 5.

## 3. Subset to the lineage of interest (optional)

If you only care about one class of cells (this protocol was developed
annotating neurons out of whole-animal data), don’t try to score every
cluster in the whole dataset against a neuron-only reference. Subset
first.

The pitfall: picking a single marker gene to threshold on. **No single
gene is expressed in every cell of a given type** — pick a gene that’s
strong in most neuron types and you will still cleanly lose whichever
types happen to express it weakly. Use a small panel of independent
markers instead, and count a cell as “in” if it clears the expression
threshold on *any* one of them (OR logic), not all of them:

``` r

# Broadly-expressed, largely-independent pan-neuronal markers -- panel,
# not singleton, specifically so that no one neuron type is missed
# because it happens to be low for any single gene.
pan_neuronal_panel <- c("unc-10", "rab-3", "ric-4", "sbt-1", "ehs-1")

expr <- Seurat::GetAssayData(obj, layer = "data")[pan_neuronal_panel, ]
is_neuronal <- Matrix::colSums(expr > 0) >= 1  # expresses at least one panel gene

neurons <- subset(obj, cells = colnames(obj)[is_neuronal])
```

Before trusting a panel like this, check it against a reference where
cell identity is already known (e.g. an existing annotated CeNGEN
dataset covering the same stage): compute the panel score for every
*known* cell type and confirm there’s a sharp break between the
lowest-scoring true-positive type and the highest-scoring true-negative
type. If the break is sharp, a simple threshold is trustworthy; if it’s
gradual, the panel needs another gene, not a fiddled threshold.

## 4. Recluster the subset

Whatever clustering resolution and PCs were used on the whole dataset
are almost never right for a subset — a lineage that was 3 clusters in
the full dataset usually needs finer resolution once it’s the only thing
in the embedding. Redo variable features, scaling, PCA, and clustering
on the subset from scratch, typically at higher `resolution`:

``` r

neurons <- Seurat::FindVariableFeatures(neurons)
neurons <- Seurat::ScaleData(neurons, features = Seurat::VariableFeatures(neurons))
neurons <- Seurat::RunPCA(neurons)
neurons <- Seurat::FindNeighbors(neurons, dims = 1:30)
neurons <- Seurat::FindClusters(neurons, resolution = 4)
neurons <- Seurat::RunUMAP(neurons, dims = 1:30)
```

## 5. Score clusters against CeNGEN

This is the part `cengenAnnotate` automates directly.

``` r

library(cengenAnnotate)

markers <- Seurat::FindAllMarkers(neurons, only.pos = TRUE)
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
