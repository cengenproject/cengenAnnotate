test_that("new_cengen_reference validates required columns", {
  expect_error(
    new_cengen_reference(tibble::tibble(gene = "g1", neuron_type = "A")),
    class = "cengen_invalid_reference_error"
  )
})

test_that("new_cengen_reference rejects duplicate (gene, neuron_type) rows", {
  bad <- tibble::tibble(
    gene = c("g1", "g1"), neuron_type = c("A", "A"),
    avg_expr = c(1, 2), pct_expr = c(10, 20)
  )
  expect_error(new_cengen_reference(bad), class = "cengen_invalid_reference_error")
})

test_that("new_cengen_reference precomputes correct per-gene mean/sd", {
  ref <- make_test_reference()
  spec_a <- ref$gene_stats[ref$gene_stats$gene == "spec_A", ]
  expect_equal(spec_a$mu, mean(c(8, 1, 1, 0, 0)), tolerance = 1e-9)
  expect_equal(spec_a$sigma, stats::sd(c(8, 1, 1, 0, 0)), tolerance = 1e-9)

  house <- ref$gene_stats[ref$gene_stats$gene == "house", ]
  expect_equal(house$sigma, 0, tolerance = 1e-9)
})

test_that("load_cengen_reference and list_cengen_datasets round-trip via a temp pins board", {
  skip_if_not_installed("pins")
  board <- pins::board_temp()

  data1 <- make_test_reference_data()
  pins::pin_write(
    board, data1,
    name = "adult_herm", type = "rds",
    title = "synthetic adult hermaphrodite",
    metadata = list(
      stage = "adult", sex = "hermaphrodite",
      dataset_label = "adult hermaphrodite (synthetic)",
      source_version = "test-fixture", date_prepared = "2026-01-01",
      n_neuron_types = length(unique(data1$neuron_type)),
      n_genes = length(unique(data1$gene))
    )
  )

  listing <- list_cengen_datasets(board = board)
  expect_equal(nrow(listing), 1)
  expect_equal(listing$name, "adult_herm")
  expect_equal(listing$stage, "adult")
  expect_equal(listing$sex, "hermaphrodite")
  expect_equal(listing$n_genes, length(unique(data1$gene)))

  ref <- load_cengen_reference("adult_herm", board = board)
  expect_s3_class(ref, "cengen_reference")
  expect_equal(nrow(ref$data), nrow(data1))
  expect_equal(ref$meta$stage, "adult")
})

test_that("list_cengen_datasets filters by pattern and handles an empty board", {
  skip_if_not_installed("pins")
  board <- pins::board_temp()

  empty_listing <- list_cengen_datasets(board = board)
  expect_equal(nrow(empty_listing), 0)

  pins::pin_write(board, make_test_reference_data(), name = "l4_herm", type = "rds")
  pins::pin_write(board, make_test_reference_data(), name = "l4_male", type = "rds")

  filtered <- list_cengen_datasets(board = board, pattern = "male")
  expect_equal(filtered$name, "l4_male")
})
