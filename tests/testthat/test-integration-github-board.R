# The only test that touches the real cengenproject/cengen-reference-data
# board over the network. Skipped by default so `devtools::test()` and CI
# stay fast and offline; run explicitly with:
#   CENGEN_RUN_INTEGRATION_TESTS=1 Rscript -e 'devtools::test()'
test_that("the live cengen-reference-data board is reachable and has at least one dataset", {
  skip_if_not(
    nzchar(Sys.getenv("CENGEN_RUN_INTEGRATION_TESTS")),
    "set CENGEN_RUN_INTEGRATION_TESTS=1 to run the live-board integration test"
  )

  board <- cengen_board()
  listing <- list_cengen_datasets(board)
  expect_gte(nrow(listing), 1)

  ref <- load_cengen_reference(listing$name[1], board = board)
  expect_s3_class(ref, "cengen_reference")
  expect_true(nrow(ref$data) > 0)
})
