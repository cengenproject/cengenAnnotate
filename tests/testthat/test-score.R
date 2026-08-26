test_that("a gene specific to one cell type scores highest there", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")

  observed <- score_cengen_matrix(panel, ref)
  expected <- expected_cluster_scores(make_test_reference_data(), "spec_A")

  observed_vec <- stats::setNames(observed$score, observed$cell_type)
  expect_equal(observed_vec[names(expected)], expected, tolerance = 1e-6)

  expect_equal(observed$cell_type[which.max(observed$score)], "A")
  expect_true(all(observed$n_genes_used == 1))
  expect_true(all(observed$n_genes_uninformative == 0))
  expect_true(all(observed$n_genes_missing_from_reference == 0))
})

test_that("a housekeeping gene (sd ~ 0) is excluded, not scored as zero, and doesn't distort other genes", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "house"))

  observed <- score_cengen_matrix(panel, ref)
  expect_true(all(observed$n_genes_uninformative == 1))
  expect_true(all(observed$n_genes_used == 1))
  expect_false(any(is.nan(observed$score)))
  expect_false(any(is.infinite(observed$score)))

  # Since "house" contributes nothing, the score must equal the single-gene
  # spec_A-only score exactly.
  solo <- score_cengen_matrix(make_panel("c1", "spec_A"), ref)
  expect_equal(observed$score, solo$score, tolerance = 1e-10)
})

test_that("coverage hinge: pct_expr just below min_pct_expr gates the score to 0", {
  ref <- make_test_reference()
  below <- score_cengen_matrix(make_panel("c1", "cov_below"), ref, min_pct_expr = 10)
  expect_equal(below$score[below$cell_type == "A"], 0, tolerance = 1e-10)
})

test_that("coverage hinge: pct_expr just above min_pct_expr yields a nonzero score", {
  ref <- make_test_reference()
  above <- score_cengen_matrix(make_panel("c1", "cov_above"), ref, min_pct_expr = 10)
  observed <- above$score[above$cell_type == "A"]
  expected <- expected_cluster_scores(make_test_reference_data(), "cov_above")[["A"]]

  expect_gt(observed, 0)
  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("a cell type with near-zero coverage across the whole panel scores exactly 0, not NaN", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "cov_below", "cov_above"))
  observed <- score_cengen_matrix(panel, ref)

  z_score <- observed$score[observed$cell_type == "Z"]
  expect_equal(z_score, 0, tolerance = 1e-10)
  expect_false(is.nan(z_score))
  expect_false(is.na(z_score))
})

test_that("a marker gene absent from the reference is dropped and tallied, other genes unaffected", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "phantom_gene"))

  observed <- score_cengen_matrix(panel, ref)
  expect_true(all(observed$n_genes_missing_from_reference == 1))
  expect_true(all(observed$n_genes_used == 1))

  solo <- score_cengen_matrix(make_panel("c1", "spec_A"), ref)
  expect_equal(observed$score, solo$score, tolerance = 1e-10)
})

test_that("a cluster whose entire panel is missing/uninformative scores NA, not an error", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("house", "phantom_gene"))

  observed <- score_cengen_matrix(panel, ref)
  expect_true(all(observed$n_genes_used == 0))
  expect_true(all(is.na(observed$score)))
})

test_that("geometric combine requires both signals; arithmetic lets one compensate the other", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "cov_above")

  geo <- score_cengen_matrix(panel, ref, combine = "geometric")
  arith <- score_cengen_matrix(panel, ref, combine = "arithmetic")

  geo_a <- geo$score[geo$cell_type == "A"]
  arith_a <- arith$score[arith$cell_type == "A"]

  # cov(A) = 0.11 (low), spec(A) is high -> geometric mean is pulled down
  # much more than the arithmetic average.
  expect_lt(geo_a, arith_a)
  expect_true(geo_a >= 0 && geo_a <= 1)
  expect_true(arith_a >= 0 && arith_a <= 1)
})

test_that("all scores are bounded in [0, 1] across a mixed panel and both combine modes", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "house", "cov_below", "cov_above", "phantom_gene"))

  for (combine in c("geometric", "arithmetic")) {
    observed <- score_cengen_matrix(panel, ref, combine = combine)
    non_na <- observed$score[!is.na(observed$score)]
    expect_true(all(non_na >= -1e-9 & non_na <= 1 + 1e-9))
  }
})

test_that("score_cengen_matrix accepts raw marker output and applies top_n/only_pos/max_p_val_adj", {
  ref <- make_test_reference()
  raw <- tibble::tibble(
    cluster = c("c1", "c1", "c1", "c1"),
    gene = c("spec_A", "cov_above", "house", "down_gene"),
    avg_log2FC = c(2.0, 1.5, 1.0, -1.0),
    p_val_adj = c(0.001, 0.001, 0.001, 0.001)
  )

  observed <- score_cengen_matrix(raw, ref, top_n = 2)
  # top_n = 2 by avg_log2FC among only_pos genes -> spec_A, cov_above
  expect_true(all(observed$n_genes_requested == 2))
})

test_that("an underscore gene symbol in the reference still matches a Seurat-renamed marker gene", {
  data <- make_test_reference_data()
  # add a gene whose "canonical" symbol has an underscore, specific to A,
  # mirroring how it would appear straight from a source table
  underscore_gene <- dplyr::mutate(
    dplyr::filter(data, gene == "spec_A"),
    gene = "CE7X_3.1"
  )
  ref <- new_cengen_reference(dplyr::bind_rows(data, underscore_gene))

  # Seurat itself would have already renamed this to a dash in the Seurat
  # object's feature names, so that's what FindAllMarkers() would report
  panel <- make_panel("c1", "CE7X-3.1")
  observed <- score_cengen_matrix(panel, ref)

  expect_true(all(observed$n_genes_missing_from_reference == 0))
  expect_equal(observed$cell_type[which.max(observed$score)], "A")
})

test_that("score_cengen_matrix scores multiple clusters independently", {
  ref <- make_test_reference()
  panel <- dplyr::bind_rows(
    make_panel("c1", "spec_A"),
    make_panel("c2", "cov_above")
  )

  observed <- score_cengen_matrix(panel, ref)
  expect_setequal(unique(observed$cluster), c("c1", "c2"))

  c1_a <- observed$score[observed$cluster == "c1" & observed$cell_type == "A"]
  c2_a <- observed$score[observed$cluster == "c2" & observed$cell_type == "A"]
  expect_false(isTRUE(all.equal(c1_a, c2_a)))
})
