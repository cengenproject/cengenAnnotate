#' Plot a CeNGEN-style dot plot for one cluster's marker panel
#'
#' Reproduces the view the CeNGEN website would show for a cluster's marker
#' gene panel, restricted to a shortlist of candidate neuron types, using
#' the same signals the coherence score is built from: dot size is percent
#' of cells expressing (coverage), dot color is the per-gene z-score
#' (specificity). Useful for visually inspecting `"review_ambiguous"`
#' clusters, or for sanity-checking an `"annotate"` call.
#'
#' @param reference A `cengen_reference` object.
#' @param markers Raw or prepared marker panel (see [score_cengen_matrix()]).
#' @param cluster The cluster to plot (matched against the `cluster`
#'   column/prepared panel).
#' @param scores Optional `cengen_scores` object (from
#'   [score_cengen_clusters()]); if supplied and `neuron_types` is not,
#'   the top `top_n_neuron_types` candidates for `cluster` are plotted.
#' @param neuron_types Optional explicit vector of neuron types to plot.
#'   Required if `scores` is not supplied.
#' @param top_n_neuron_types Number of top-scoring neuron types to plot
#'   when deriving the shortlist from `scores` (default 8).
#' @param top_n,cluster_col,gene_col Passed to [prepare_marker_panel()] when
#'   `markers` is not already a prepared panel.
#'
#' @return A `ggplot` object.
#' @export
plot_cengen_dotplot <- function(
  reference,
  markers,
  cluster,
  scores = NULL,
  neuron_types = NULL,
  top_n_neuron_types = 8,
  top_n = 20,
  cluster_col = "cluster",
  gene_col = "gene"
) {
  panel <- if (is_prepared_panel(markers)) {
    tibble::as_tibble(markers)
  } else {
    prepare_marker_panel(markers, top_n = top_n, cluster_col = cluster_col, gene_col = gene_col)
  }
  panel$gene <- normalize_gene_symbol(as.character(panel$gene))

  cluster <- as.character(cluster)
  cl_panel <- panel[panel$cluster == cluster, , drop = FALSE]
  if (nrow(cl_panel) == 0) {
    cengen_abort(sprintf("No marker genes found for cluster '%s'.", cluster), class = "cengen_no_markers_error")
  }
  genes <- unique(cl_panel$gene)

  if (is.null(neuron_types)) {
    if (is.null(scores)) {
      cengen_abort(
        "Provide either `neuron_types` or `scores` (from score_cengen_clusters()) to choose which neuron types to plot.",
        class = "cengen_missing_neuron_types_error"
      )
    }
    mat <- cengen_matrix(scores)
    if (is.null(mat)) {
      cengen_abort(
        "`scores` has no attached score matrix; pass the object returned by score_cengen_clusters(), or supply `neuron_types` directly.",
        class = "cengen_missing_score_matrix_error"
      )
    }
    sub <- mat[mat$cluster == cluster, , drop = FALSE]
    sub <- sub[order(-ifelse(is.na(sub$score), -Inf, sub$score)), , drop = FALSE]
    neuron_types <- utils::head(sub$neuron_type, top_n_neuron_types)
  }

  plot_data <- reference$data[
    reference$data$gene %in% genes & reference$data$neuron_type %in% neuron_types,
    , drop = FALSE
  ]
  plot_data <- dplyr::left_join(plot_data, reference$gene_stats, by = "gene")
  plot_data$z <- (plot_data$avg_expr - plot_data$mu) / plot_data$sigma
  plot_data$gene <- factor(plot_data$gene, levels = rev(genes))
  plot_data$neuron_type <- factor(plot_data$neuron_type, levels = neuron_types)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$neuron_type, y = .data$gene, size = .data$pct_expr, color = .data$z)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "grey85", high = "#B2182B", name = "z-score") +
    ggplot2::scale_size_continuous(name = "% expressing", range = c(0, 8), limits = c(0, 100)) +
    ggplot2::labs(x = "Neuron type", y = "Marker gene", title = paste("Cluster", cluster)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
