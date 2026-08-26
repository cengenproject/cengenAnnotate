make_test_seurat <- function(n_genes = 5, n_cells = 20, n_clusters = 2) {
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("cell", seq_len(n_cells)))
  )
  obj <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  )
  obj$seurat_clusters <- factor(rep(
    as.character(seq_len(n_clusters) - 1L),
    length.out = n_cells
  ))
  obj
}
