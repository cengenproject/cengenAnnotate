test_that("score_cengen_clusters wraps scoring + classification and attaches the full matrix", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")

  result <- score_cengen_clusters(panel, ref)
  expect_s3_class(result, "cengen_scores")
  expect_equal(result$cluster, "c1")
  expect_equal(result$best_neuron_type, "A")

  matrix <- cengen_matrix(result)
  expect_true(!is.null(matrix))
  expect_setequal(matrix$neuron_type, c("A", "B", "C", "D", "Z"))

  direct <- classify_cengen_matches(score_cengen_matrix(panel, ref))
  expect_equal(result$decision, direct$decision)
})

test_that("cengen_matrix() returns NULL for an object with no attached matrix", {
  expect_null(cengen_matrix(tibble::tibble(cluster = "c1")))
})
