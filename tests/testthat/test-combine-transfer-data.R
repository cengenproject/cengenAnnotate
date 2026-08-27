make_td_test_object <- function() {
  # 4 clusters x 5 cells:
  #   "0": cengenAnnotate ambiguous, TransferData confident -> should become provisional
  #   "1": cengenAnnotate ambiguous, TransferData NOT confident -> stays ambiguous
  #   "2": cengenAnnotate already annotate, TransferData agrees -> stays annotate
  #   "3": cengenAnnotate already annotate, TransferData confidently disagrees -> conflict
  obj <- make_test_seurat(n_genes = 5, n_cells = 20, n_clusters = 4)
  cl <- as.character(obj$seurat_clusters)

  predicted <- character(length(cl))
  score <- numeric(length(cl))

  idx0 <- which(cl == "0"); predicted[idx0] <- c("AWA", "AWA", "AWA", "AWA", "ASH"); score[idx0] <- 0.9
  idx1 <- which(cl == "1"); predicted[idx1] <- c("AIY", "RIA", "URX", "AWB", "AWA"); score[idx1] <- 0.9
  idx2 <- which(cl == "2"); predicted[idx2] <- rep("FOO", 5); score[idx2] <- 0.9
  idx3 <- which(cl == "3"); predicted[idx3] <- rep("BAR", 5); score[idx3] <- 0.9

  obj$predicted.id <- predicted
  obj$prediction.score.max <- score
  obj
}

make_td_test_scores <- function() {
  tibble::tribble(
    ~cluster, ~best_cell_type, ~best_score, ~decision,
    "0", "whatever", 0.4, "review_ambiguous",
    "1", "whatever", 0.4, "review_ambiguous",
    "2", "FOO", 0.9, "annotate",
    "3", "BAZ", 0.9, "annotate"
  )
}

test_that("a TransferData-confident cluster cengenAnnotate left ambiguous becomes provisional", {
  scores <- make_td_test_scores()
  obj <- make_td_test_object()
  result <- combine_with_transfer_data(scores, obj)

  row0 <- result[result$cluster == "0", ]
  expect_equal(row0$decision, "provisional")
  expect_equal(row0$best_cell_type, "AWA")
  expect_equal(row0$td_call, "AWA")
  expect_equal(row0$td_purity, 0.8)
})

test_that("a cluster where TransferData is not confident keeps its original decision", {
  scores <- make_td_test_scores()
  obj <- make_td_test_object()
  result <- combine_with_transfer_data(scores, obj)

  row1 <- result[result$cluster == "1", ]
  expect_equal(row1$decision, "review_ambiguous")
  expect_equal(row1$best_cell_type, "whatever")
})

test_that("an already-annotated cluster where TransferData agrees stays annotate", {
  scores <- make_td_test_scores()
  obj <- make_td_test_object()
  result <- combine_with_transfer_data(scores, obj)

  row2 <- result[result$cluster == "2", ]
  expect_equal(row2$decision, "annotate")
  expect_equal(row2$best_cell_type, "FOO")
  expect_equal(row2$td_call, "FOO")
})

test_that("an already-annotated cluster where TransferData confidently disagrees is flagged review_conflict", {
  scores <- make_td_test_scores()
  obj <- make_td_test_object()
  result <- combine_with_transfer_data(scores, obj)

  row3 <- result[result$cluster == "3", ]
  expect_equal(row3$decision, "review_conflict")
  # cengenAnnotate's own (independently validated) call is preserved, not overwritten
  expect_equal(row3$best_cell_type, "BAZ")
  expect_equal(row3$td_call, "BAR")
})

test_that("thresholds are respected - raising td_min_purity un-confidents cluster 0", {
  scores <- make_td_test_scores()
  obj <- make_td_test_object()
  result <- combine_with_transfer_data(scores, obj, td_min_purity = 0.9)

  row0 <- result[result$cluster == "0", ]
  expect_equal(row0$decision, "review_ambiguous")
  expect_equal(row0$best_cell_type, "whatever")
})

test_that("the score matrix attribute is preserved through combine_with_transfer_data()", {
  scores <- make_td_test_scores()
  attr(scores, "matrix") <- tibble::tibble(cluster = "0", cell_type = "AWA", score = 0.9)
  class(scores) <- c("cengen_scores", class(scores))
  obj <- make_td_test_object()

  result <- combine_with_transfer_data(scores, obj)
  expect_s3_class(result, "cengen_scores")
  expect_false(is.null(cengen_matrix(result)))
  expect_equal(cengen_matrix(result)$cell_type, "AWA")
})

test_that("combine_with_transfer_data requires a Seurat object", {
  expect_error(
    combine_with_transfer_data(make_td_test_scores(), data.frame(x = 1)),
    class = "cengen_invalid_object_error"
  )
})

test_that("combine_with_transfer_data errors informatively on missing TransferData columns", {
  obj <- make_test_seurat(n_genes = 5, n_cells = 10, n_clusters = 2)
  expect_error(
    combine_with_transfer_data(make_td_test_scores(), obj),
    class = "cengen_invalid_object_error"
  )
})
