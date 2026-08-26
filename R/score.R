#' Score clusters against every candidate CeNGEN cell type
#'
#' Computes a "dot-plot-mirrored coherence score" for each (cluster,
#' cell type) pair: how consistently the cluster's marker gene panel is
#' both highly expressed (coverage, mirroring dot *size*) and specific
#' (mirroring dot *color*) in that cell type, relative to CeNGEN's other
#' cell types. This is the computational equivalent of visually judging
#' whether a CeNGEN dot plot shows one cell type standing out.
#'
#' For each gene `g` and cell type `t`:
#' \itemize{
#'   \item `cov(g,t) = 0` if `pct_expr(g,t) < min_pct_expr`, else
#'     `pct_expr(g,t) / 100`.
#'   \item `z(g,t) = (avg_expr(g,t) - mu_g) / sigma_g`, where `mu_g`/`sigma_g`
#'     are gene `g`'s mean/sd of `avg_expr` across all cell types.
#'     `spec(g,t) = plogis(z(g,t))`.
#'   \item `s(g,t)` combines the two: for `combine = "geometric"` (default),
#'     `cov(g,t)^cov_weight * spec(g,t)^spec_weight` (requires both signals
#'     to be non-trivial simultaneously); for `"arithmetic"`,
#'     `cov_weight * cov(g,t) + spec_weight * spec(g,t)`.
#' }
#' The cluster-level score for cell type `t` is the mean of `s(g,t)` over
#' the cluster's usable marker genes. Genes absent from the reference, or
#' with near-zero variance across cell types (uninformative), are
#' excluded from scoring and tallied in the diagnostic columns.
#'
#' @param markers Either raw marker output (see [prepare_marker_panel()])
#'   or an already-prepared panel (a data frame with `cluster`, `gene`,
#'   `rank`, `n_used` columns).
#' @param reference A `cengen_reference` object (see [new_cengen_reference()]
#'   / [load_cengen_reference()]).
#' @param top_n,cluster_col,gene_col Passed to [prepare_marker_panel()] when
#'   `markers` is not already a prepared panel.
#' @param min_pct_expr Minimum percent-expressing (0-100) for a gene to
#'   count as "covered" in a cell type; below this, coverage is 0.
#' @param combine How to combine the coverage and specificity signals per
#'   gene: `"geometric"` (default) or `"arithmetic"`.
#' @param cov_weight,spec_weight Weights on the coverage and specificity
#'   signals (exponents for `"geometric"`, linear weights for
#'   `"arithmetic"`).
#'
#' @return A long tibble with columns `cluster`, `cell_type`, `score`,
#'   `n_genes_used`, `n_genes_requested`, `n_genes_missing_from_reference`,
#'   `n_genes_uninformative`.
#' @export
score_cengen_matrix <- function(
  markers,
  reference,
  top_n = 20,
  min_pct_expr = 10,
  combine = c("geometric", "arithmetic"),
  cov_weight = 0.5,
  spec_weight = 0.5,
  cluster_col = "cluster",
  gene_col = "gene"
) {
  combine <- match.arg(combine)
  if (!inherits(reference, "cengen_reference")) {
    cengen_abort(
      "`reference` must be a cengen_reference object (see new_cengen_reference() / load_cengen_reference()).",
      class = "cengen_invalid_reference_error"
    )
  }

  panel <- if (is_prepared_panel(markers)) {
    tibble::as_tibble(markers)
  } else {
    prepare_marker_panel(markers, top_n = top_n, cluster_col = cluster_col, gene_col = gene_col)
  }
  # Normalize unconditionally (idempotent if already normalized by
  # prepare_marker_panel()) so a pre-built panel handed in directly is held
  # to the same convention as the reference's gene column.
  panel$gene <- normalize_gene_symbol(as.character(panel$gene))

  cell_types <- sort(unique(reference$data$cell_type))
  ref_genes <- unique(reference$gene_stats$gene)
  eps <- 1e-8

  clusters <- unique(panel$cluster)

  results <- lapply(clusters, function(cl) {
    cl_panel <- panel[panel$cluster == cl, , drop = FALSE]
    genes <- unique(cl_panel$gene)
    n_requested <- cl_panel$n_used[1]

    missing_genes <- setdiff(genes, ref_genes)
    present_genes <- intersect(genes, ref_genes)

    gs <- reference$gene_stats[reference$gene_stats$gene %in% present_genes, , drop = FALSE]
    uninformative_genes <- gs$gene[is.na(gs$sigma) | gs$sigma <= eps]
    usable_genes <- setdiff(present_genes, uninformative_genes)

    n_missing <- length(missing_genes)
    n_uninformative <- length(uninformative_genes)
    n_used <- length(usable_genes)

    if (n_used == 0) {
      return(tibble::tibble(
        cluster = cl,
        cell_type = cell_types,
        score = NA_real_,
        n_genes_used = n_used,
        n_genes_requested = n_requested,
        n_genes_missing_from_reference = n_missing,
        n_genes_uninformative = n_uninformative
      ))
    }

    sub <- reference$data[reference$data$gene %in% usable_genes, , drop = FALSE]
    gs_sub <- reference$gene_stats[reference$gene_stats$gene %in% usable_genes, , drop = FALSE]
    sub <- dplyr::left_join(sub, gs_sub, by = "gene")

    sub$cov <- ifelse(sub$pct_expr < min_pct_expr, 0, sub$pct_expr / 100)
    sub$z <- (sub$avg_expr - sub$mu) / sub$sigma
    sub$spec <- stats::plogis(sub$z)

    sub$s <- if (combine == "geometric") {
      (sub$cov^cov_weight) * (sub$spec^spec_weight)
    } else {
      cov_weight * sub$cov + spec_weight * sub$spec
    }

    agg <- dplyr::summarise(
      dplyr::group_by(sub, .data$cell_type),
      score = mean(.data$s, na.rm = TRUE),
      .groups = "drop"
    )

    full <- tibble::tibble(cell_type = cell_types)
    agg <- dplyr::left_join(full, agg, by = "cell_type")

    tibble::tibble(
      cluster = cl,
      cell_type = agg$cell_type,
      score = agg$score,
      n_genes_used = n_used,
      n_genes_requested = n_requested,
      n_genes_missing_from_reference = n_missing,
      n_genes_uninformative = n_uninformative
    )
  })

  dplyr::bind_rows(results)
}
