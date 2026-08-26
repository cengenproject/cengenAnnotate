test_that("plot_cengen_dotplot returns a ggplot when cell_types is given explicitly", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "cov_above"))

  p <- plot_cengen_dotplot(ref, panel, "c1", cell_types = c("A", "B", "Z"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_cengen_dotplot derives the cell type shortlist from scores", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "cov_above"))
  scores <- score_cengen_clusters(panel, ref)

  p <- plot_cengen_dotplot(ref, panel, "c1", scores = scores, top_n_cell_types = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cengen_dotplot errors without cell_types or scores", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")
  expect_error(
    plot_cengen_dotplot(ref, panel, "c1"),
    class = "cengen_missing_cell_types_error"
  )
})

test_that("plot_cengen_dotplot errors for an unknown cluster", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")
  expect_error(
    plot_cengen_dotplot(ref, panel, "does_not_exist", cell_types = "A"),
    class = "cengen_no_markers_error"
  )
})
