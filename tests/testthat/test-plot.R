test_that("plot_cengen_dotplot returns a ggplot when neuron_types is given explicitly", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "cov_above"))

  p <- plot_cengen_dotplot(ref, panel, "c1", neuron_types = c("A", "B", "Z"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_cengen_dotplot derives the neuron type shortlist from scores", {
  ref <- make_test_reference()
  panel <- make_panel("c1", c("spec_A", "cov_above"))
  scores <- score_cengen_clusters(panel, ref)

  p <- plot_cengen_dotplot(ref, panel, "c1", scores = scores, top_n_neuron_types = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cengen_dotplot errors without neuron_types or scores", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")
  expect_error(
    plot_cengen_dotplot(ref, panel, "c1"),
    class = "cengen_missing_neuron_types_error"
  )
})

test_that("plot_cengen_dotplot errors for an unknown cluster", {
  ref <- make_test_reference()
  panel <- make_panel("c1", "spec_A")
  expect_error(
    plot_cengen_dotplot(ref, panel, "does_not_exist", neuron_types = "A"),
    class = "cengen_no_markers_error"
  )
})
