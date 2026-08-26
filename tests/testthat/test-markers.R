make_raw_markers <- function() {
  tibble::tribble(
    ~cluster, ~gene, ~avg_log2FC, ~p_val_adj,
    "c1", "g1", 3.0, 0.001,
    "c1", "g2", 2.5, 0.001,
    "c1", "g3", 2.0, 0.001,
    "c1", "g4", 1.5, 0.001,
    "c1", "g5", 1.0, 0.001,
    "c2", "h1", 4.0, 0.001,
    "c2", "h2", 3.0, 0.20, # filtered out by p_val_adj
    "c2", "h3", -1.0, 0.001 # filtered out by only_pos
  )
}

test_that("more markers than top_n are truncated to the top N by fold-change", {
  panel <- prepare_marker_panel(make_raw_markers(), top_n = 3)
  c1 <- panel[panel$cluster == "c1", ]
  expect_equal(nrow(c1), 3)
  expect_equal(c1$gene, c("g1", "g2", "g3"))
  expect_equal(c1$rank, 1:3)
  expect_true(all(c1$n_used == 3))
  expect_true(all(c1$n_requested == 3))
})

test_that("exactly top_n available markers uses all of them", {
  panel <- prepare_marker_panel(make_raw_markers(), top_n = 5)
  c1 <- panel[panel$cluster == "c1", ]
  expect_equal(nrow(c1), 5)
  expect_true(all(c1$n_used == 5))
})

test_that("fewer than top_n available markers uses only what's available, and reports n_used < top_n", {
  panel <- prepare_marker_panel(make_raw_markers(), top_n = 20)
  c2 <- panel[panel$cluster == "c2", ]
  # only h1 survives only_pos + max_p_val_adj filtering
  expect_equal(nrow(c2), 1)
  expect_equal(c2$gene, "h1")
  expect_true(all(c2$n_used == 1))
  expect_true(all(c2$n_requested == 20))
})

test_that("only_pos = FALSE keeps negative fold-change genes", {
  panel <- prepare_marker_panel(make_raw_markers(), top_n = 20, only_pos = FALSE, max_p_val_adj = NULL)
  c2 <- panel[panel$cluster == "c2", ]
  expect_true("h3" %in% c2$gene)
})

test_that("max_p_val_adj = NULL skips the significance filter", {
  panel <- prepare_marker_panel(make_raw_markers(), top_n = 20, max_p_val_adj = NULL)
  c2 <- panel[panel$cluster == "c2", ]
  expect_true("h2" %in% c2$gene)
})

test_that("duplicate gene rows within a cluster are deduplicated", {
  raw <- tibble::tribble(
    ~cluster, ~gene, ~avg_log2FC, ~p_val_adj,
    "c1", "g1", 3.0, 0.001,
    "c1", "g1", 2.0, 0.001,
    "c1", "g2", 1.0, 0.001
  )
  panel <- prepare_marker_panel(raw, top_n = 20)
  expect_equal(sum(panel$gene == "g1"), 1)
})

test_that("gene symbols with underscores are normalized to dashes (Seurat's own feature-name convention)", {
  raw <- tibble::tribble(
    ~cluster, ~gene, ~avg_log2FC, ~p_val_adj,
    "c1", "CE7X_3.1", 3.0, 0.001
  )
  panel <- prepare_marker_panel(raw, top_n = 20)
  expect_equal(panel$gene, "CE7X-3.1")
})

test_that("missing required columns error informatively", {
  expect_error(
    prepare_marker_panel(tibble::tibble(cluster = "c1", gene = "g1")),
    class = "cengen_invalid_markers_error"
  )
})
