make_score_row <- function(cluster, cell_type, score, n_genes_used = 10,
                            n_genes_requested = 10, n_genes_missing = 0, n_genes_uninformative = 0) {
  tibble::tibble(
    cluster = cluster, cell_type = cell_type, score = score,
    n_genes_used = n_genes_used, n_genes_requested = n_genes_requested,
    n_genes_missing_from_reference = n_genes_missing,
    n_genes_uninformative = n_genes_uninformative
  )
}

test_that("a clear winner (high score, large gap) is annotated", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.80),
    make_score_row("c1", "B", 0.30),
    make_score_row("c1", "C", 0.10)
  )
  result <- classify_cengen_matches(m)
  expect_equal(result$decision, "annotate")
  expect_equal(result$best_cell_type, "A")
  expect_equal(result$second_cell_type, "B")
  expect_equal(result$gap, 0.50, tolerance = 1e-9)
})

test_that("a high score with a small gap to second place is ambiguous", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.70),
    make_score_row("c1", "B", 0.65),
    make_score_row("c1", "C", 0.10)
  )
  result <- classify_cengen_matches(m, cutoff = 0.6, min_gap = 0.10)
  expect_equal(result$decision, "review_ambiguous")
})

test_that("a low top score is ambiguous even with a large gap", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.40),
    make_score_row("c1", "B", 0.05)
  )
  result <- classify_cengen_matches(m, cutoff = 0.6, min_gap = 0.10)
  expect_equal(result$decision, "review_ambiguous")
})

test_that("insufficient marker genes overrides an otherwise-clear match", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.90, n_genes_used = 2),
    make_score_row("c1", "B", 0.10, n_genes_used = 2)
  )
  result <- classify_cengen_matches(m, min_genes_used = 5)
  expect_equal(result$decision, "review_insufficient_markers")
})

test_that("boundary values (score == cutoff, gap == min_gap) are inclusive and annotate", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.60),
    make_score_row("c1", "B", 0.50)
  )
  result <- classify_cengen_matches(m, cutoff = 0.6, min_gap = 0.10)
  expect_equal(result$decision, "annotate")
})

test_that("a single-candidate matrix (no second cell type) never auto-annotates", {
  m <- make_score_row("c1", "A", 0.95)
  result <- classify_cengen_matches(m, cutoff = 0.6, min_gap = 0.10)
  expect_true(is.na(result$second_cell_type))
  expect_true(is.na(result$gap))
  expect_equal(result$decision, "review_ambiguous")
})

test_that("classify_cengen_matches handles multiple clusters independently", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", 0.80), make_score_row("c1", "B", 0.10),
    make_score_row("c2", "A", 0.30), make_score_row("c2", "B", 0.20)
  )
  result <- classify_cengen_matches(m)
  expect_setequal(result$cluster, c("c1", "c2"))
  expect_equal(result$decision[result$cluster == "c1"], "annotate")
  expect_equal(result$decision[result$cluster == "c2"], "review_ambiguous")
})

test_that("results are sorted by best_score descending, not by input cluster order", {
  m <- dplyr::bind_rows(
    make_score_row("low_conf", "A", 0.20), make_score_row("low_conf", "B", 0.05),
    make_score_row("high_conf", "A", 0.90), make_score_row("high_conf", "B", 0.10),
    make_score_row("mid_conf", "A", 0.55), make_score_row("mid_conf", "B", 0.20)
  )
  result <- classify_cengen_matches(m)
  expect_equal(result$cluster, c("high_conf", "mid_conf", "low_conf"))
})

test_that("clusters with an NA best_score (all markers excluded) sort last, not erroring", {
  m <- dplyr::bind_rows(
    make_score_row("c1", "A", NA_real_, n_genes_used = 0),
    make_score_row("c2", "A", 0.5)
  )
  result <- classify_cengen_matches(m, min_genes_used = 0)
  expect_equal(result$cluster, c("c2", "c1"))
})

test_that("classify_cengen_matches errors informatively on a malformed score matrix", {
  expect_error(
    classify_cengen_matches(tibble::tibble(cluster = "c1", cell_type = "A")),
    class = "cengen_invalid_score_matrix_error"
  )
})
