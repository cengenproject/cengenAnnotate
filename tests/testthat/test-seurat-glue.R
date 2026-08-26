make_test_scores <- function() {
  tibble::tribble(
    ~cluster, ~best_neuron_type, ~best_score, ~decision,
    "0", "AWA", 0.85, "annotate",
    "1", "AVA", 0.40, "review_ambiguous"
  )
}

test_that("write_cengen_annotations writes only 'annotate' clusters by default", {
  obj <- make_test_seurat()
  result <- suppressWarnings(write_cengen_annotations(obj, make_test_scores()))

  ann <- result@meta.data$cengen_annotation
  clusters <- as.character(result@meta.data$seurat_clusters)

  expect_true(all(ann[clusters == "0"] == "AWA"))
  expect_true(all(is.na(ann[clusters == "1"])))
})

test_that("write_cengen_annotations also writes the matched score", {
  obj <- make_test_seurat()
  result <- suppressWarnings(write_cengen_annotations(obj, make_test_scores()))

  scr <- result@meta.data$cengen_annotation_score
  clusters <- as.character(result@meta.data$seurat_clusters)
  expect_true(all(scr[clusters == "0"] == 0.85))
  expect_true(all(is.na(scr[clusters == "1"])))
})

test_that("include widens which decisions get written", {
  obj <- make_test_seurat()
  result <- suppressWarnings(write_cengen_annotations(
    obj, make_test_scores(),
    include = c("annotate", "review_ambiguous")
  ))
  ann <- result@meta.data$cengen_annotation
  clusters <- as.character(result@meta.data$seurat_clusters)
  expect_true(all(ann[clusters == "1"] == "AVA"))
})

test_that("overwrite = FALSE (default) errors if the metadata column already exists", {
  obj <- make_test_seurat()
  once <- suppressWarnings(write_cengen_annotations(obj, make_test_scores()))
  expect_error(
    write_cengen_annotations(once, make_test_scores()),
    class = "cengen_metadata_exists_error"
  )
  expect_no_error(
    suppressWarnings(write_cengen_annotations(once, make_test_scores(), overwrite = TRUE))
  )
})

test_that("a cluster in `scores` with no match in `object` warns", {
  obj <- make_test_seurat()
  scores <- make_test_scores()
  scores$cluster[2] <- "not_a_real_cluster"
  expect_warning(
    expect_warning(write_cengen_annotations(obj, scores), regexp = "no matching row in `scores`"),
    regexp = "no matching cluster in `object`"
  )
})

test_that("a cluster in `object` with no match in `scores` warns", {
  obj <- make_test_seurat()
  scores <- make_test_scores()[1, ]
  expect_warning(write_cengen_annotations(obj, scores), regexp = "no matching row in `scores`")
})

test_that("write_cengen_annotations requires a Seurat object", {
  expect_error(
    write_cengen_annotations(data.frame(x = 1), make_test_scores()),
    class = "cengen_invalid_object_error"
  )
})
