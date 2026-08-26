# Small synthetic CeNGEN-like reference used across scoring tests.
#
# Cell types: A, B, C, D, Z
#   - "spec_A" is specific to A (high avg_expr + high pct_expr there, low
#     elsewhere), and reasonably covered (pct_expr >= 10) everywhere except Z.
#   - "house" is a housekeeping gene: identical avg_expr everywhere, so its
#     across-type sd is 0 -> must be excluded from scoring (uninformative).
#   - "cov_below"/"cov_above" share the same avg_expr pattern (specific to A)
#     but differ only in pct_expr at A: 9 (just below the default
#     min_pct_expr = 10) vs 11 (just above), isolating the coverage hinge.
#   - Cell type Z has pct_expr = 2 (< 10) for every gene, so every gene's
#     coverage is gated to 0 there -> Z's score is always exactly 0, with no
#     dependence on specificity.
make_test_reference_data <- function() {
  tibble::tribble(
    ~gene,        ~cell_type, ~avg_expr, ~pct_expr,
    "spec_A",     "A",          8,         95,
    "spec_A",     "B",          1,         20,
    "spec_A",     "C",          1,         20,
    "spec_A",     "D",          0,         20,
    "spec_A",     "Z",          0,         2,

    "house",      "A",          5,         80,
    "house",      "B",          5,         80,
    "house",      "C",          5,         80,
    "house",      "D",          5,         80,
    "house",      "Z",          5,         2,

    "cov_below",  "A",          6,         9,
    "cov_below",  "B",          0,         5,
    "cov_below",  "C",          0,         5,
    "cov_below",  "D",          0,         5,
    "cov_below",  "Z",          0,         2,

    "cov_above",  "A",          6,         11,
    "cov_above",  "B",          0,         5,
    "cov_above",  "C",          0,         5,
    "cov_above",  "D",          0,         5,
    "cov_above",  "Z",          0,         2
  )
}

make_panel <- function(cluster, genes) {
  tibble::tibble(
    cluster = cluster, gene = genes, rank = seq_along(genes),
    n_requested = length(genes), n_used = length(genes)
  )
}

make_test_reference <- function() {
  new_cengen_reference(
    make_test_reference_data(),
    meta = list(stage = "test", sex = "test", dataset_label = "synthetic test reference")
  )
}

# Independent (non-package) reference implementation of the per-gene,
# per-cell-type coherence signal, used to compute expected values without
# sharing code paths with R/score.R.
expected_gene_type_signal <- function(data, gene, min_pct_expr = 10,
                                       combine = "geometric",
                                       cov_weight = 0.5, spec_weight = 0.5) {
  rows <- data[data$gene == gene, , drop = FALSE]
  mu <- mean(rows$avg_expr)
  sigma <- stats::sd(rows$avg_expr)

  cov <- ifelse(rows$pct_expr < min_pct_expr, 0, rows$pct_expr / 100)
  z <- (rows$avg_expr - mu) / sigma
  spec <- stats::plogis(z)

  s <- if (combine == "geometric") {
    (cov^cov_weight) * (spec^spec_weight)
  } else {
    cov_weight * cov + spec_weight * spec
  }

  stats::setNames(s, rows$cell_type)
}

# Expected cluster-level score per cell type for a given gene set, using
# the same exclusion rules as the package (missing genes and sigma ~ 0
# genes are dropped before averaging).
expected_cluster_scores <- function(data, genes, min_pct_expr = 10,
                                     combine = "geometric",
                                     cov_weight = 0.5, spec_weight = 0.5,
                                     eps = 1e-8) {
  present <- intersect(genes, unique(data$gene))
  usable <- Filter(function(g) {
    rows <- data[data$gene == g, , drop = FALSE]
    sigma <- stats::sd(rows$avg_expr)
    !is.na(sigma) && sigma > eps
  }, present)

  cell_types <- sort(unique(data$cell_type))
  if (length(usable) == 0) {
    return(stats::setNames(rep(NA_real_, length(cell_types)), cell_types))
  }

  signals <- lapply(usable, expected_gene_type_signal,
    data = data, min_pct_expr = min_pct_expr, combine = combine,
    cov_weight = cov_weight, spec_weight = spec_weight
  )
  mat <- do.call(rbind, signals)
  colMeans(mat)[cell_types]
}
