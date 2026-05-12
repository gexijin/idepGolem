# Unit tests for fct_03_clustering.R
#
# Covers the bug fixes made for GitHub issue #856:
#   "Cluster IDs are messed up in k-means clustering enrichment analysis"
#
# Bugs fixed:
#   1. set.seed() moved to immediately before draw() — reproducibility
#   2. kmeans_gene_lists used shiny_env$ht@matrix (no slot) — fixed by storing
#      shiny_env$ht_matrix separately at render time
#   3. current_config omitted heatmap_cutoff/gene_centering/gene_normalize —
#      causing stale cluster IDs after those settings changed
#   4. heatmap_main_object crashed on transient n_genes = 0 while typing
#
# Functions tested directly (no Shiny reactive context needed):
#   process_heatmap_data()
#   heatmap_main()        — k-means cluster ID names and reproducibility


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Synthetic expression matrix: 150 genes x 20 samples
set.seed(1)
test_mat <- matrix(rnorm(150 * 20), nrow = 150, ncol = 20)
rownames(test_mat) <- paste0("GENE", formatC(seq_len(150), width = 3, flag = "0"))
colnames(test_mat) <- paste0("S", seq_len(20))

# Minimal all_gene_names: one column triggers the "User_ID" passthrough in
# rowname_id_swap() so row names are left untouched.
dummy_gene_names <- data.frame(
  dummy = seq_len(nrow(test_mat)),
  row.names = rownames(test_mat)
)

# Pre-processed matrix used as input to heatmap_main()
hmap_data <- process_heatmap_data(
  data             = test_mat,
  n_genes_max      = 100,
  gene_centering   = TRUE,
  gene_normalize   = FALSE,
  sample_centering = FALSE,
  sample_normalize = FALSE,
  all_gene_names   = dummy_gene_names,
  select_gene_id   = "User_ID"
)

# Helper: call heatmap_main() with k-means defaults, drawing to a temp PDF
# so no interactive graphics device is required.
# Automatically skips the calling test if ComplexHeatmap is not installed.
run_kmeans_heatmap <- function(k, seed = 0, data = hmap_data) {
  skip_if_not_installed("ComplexHeatmap")
  dist_funs   <- dist_functions()
  hclust_funs <- hcluster_functions()

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)

  heatmap_main(
    data                   = data,
    cluster_meth           = 2,           # k-means
    heatmap_cutoff         = 3,
    sample_info            = NULL,
    select_factors_heatmap = "None",
    dist_funs              = dist_funs,
    dist_function          = "1",         # Pearson (unused for k-means)
    hclust_function        = "average",   # unused for k-means
    sample_clustering      = FALSE,
    heatmap_color_select   = c("blue", "white", "red"),
    row_dend               = FALSE,
    k_clusters             = k,
    re_run                 = seed,
    selected_genes         = NULL,        # NULL avoids req() in non-reactive ctx
    group_pal              = NULL,
    sample_color           = "Okabe Ito",
    show_column_names      = TRUE,
    show_heatmap_legend    = FALSE,
    row_dend_obj           = NULL,
    col_dend_obj           = NULL,
    use_letter_overlay     = FALSE,
    show_cluster_labels    = FALSE,
    custom_cluster_labels  = NULL
  )
}


# ===========================================================================
# SECTION A — process_heatmap_data() edge cases
# ===========================================================================

test_that("A1: process_heatmap_data returns n_genes rows with gene_centering=TRUE", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 80,
    gene_centering = TRUE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  expect_equal(nrow(result), 80)
  expect_equal(ncol(result), ncol(test_mat))
})

test_that("A2: process_heatmap_data clamps n_genes to nrow(data) when too large", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 99999,
    gene_centering = TRUE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  expect_lte(nrow(result), nrow(test_mat))
})

test_that("A3: process_heatmap_data clamps n_genes to 10 when below minimum", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 3,        # below minimum of 10
    gene_centering = TRUE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  expect_equal(nrow(result), 10)
})

test_that("A4: process_heatmap_data with gene_centering=FALSE still returns n_genes rows", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 50,
    gene_centering = FALSE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  expect_equal(nrow(result), 50)
})

test_that("A5: process_heatmap_data with gene_normalize=TRUE produces near-unit row SDs", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 50,
    gene_centering = TRUE,
    gene_normalize = TRUE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  row_sds <- apply(result, 1, sd)
  # All row SDs should be close to 1 after normalization
  # (small tolerance for rounding)
  expect_true(all(abs(row_sds - 1) < 0.01),
    info = paste("Row SDs:", paste(round(row_sds[1:5], 3), collapse = ", ")))
})


# ===========================================================================
# SECTION B — k-means cluster ID names (core bug #856 regression)
#
# row_order() on a drawn k-means HeatmapList returns a named list where
# NAMES are the cluster numbers as character strings "1", "2", ..., "k".
# Before the fix, these could be inflated row-index values (e.g., "500",
# "2000") because kmeans_gene_lists was reading from a mismatched matrix.
# ===========================================================================

test_that("B1: k-means row_order names are exactly '1' through k for k=5", {
  ht <- run_kmeans_heatmap(k = 5, seed = 42)

  expect_s4_class(ht, "HeatmapList")

  ro <- ComplexHeatmap::row_order(ht)
  expect_true(is.list(ro), info = "row_order() should return a list for k-means")
  expect_equal(length(ro), 5)

  # Names must be the cluster numbers "1":"5", not inflated row indices
  expect_equal(sort(as.integer(names(ro))), 1:5,
    info = paste("Actual names:", paste(names(ro), collapse = ", ")))
})

test_that("B2: k-means row_order names are exactly '1' through '10' for k=10", {
  ht <- run_kmeans_heatmap(k = 10, seed = 7)

  ro <- ComplexHeatmap::row_order(ht)
  expect_equal(length(ro), 10)
  expect_equal(sort(as.integer(names(ro))), 1:10,
    info = paste("Actual names:", paste(names(ro), collapse = ", ")))
})

test_that("B3: every row index in row_order is within valid matrix row range", {
  k   <- 6
  ht  <- run_kmeans_heatmap(k = k, seed = 0)
  ro  <- ComplexHeatmap::row_order(ht)
  n   <- nrow(hmap_data)

  all_indices <- unlist(ro)
  expect_true(
    all(all_indices >= 1 & all_indices <= n),
    info = sprintf(
      "Row indices out of range [1, %d]: %s",
      n, paste(all_indices[all_indices < 1 | all_indices > n], collapse = ", ")
    )
  )
})

test_that("B4: all k clusters are non-empty (no empty slices)", {
  ht <- run_kmeans_heatmap(k = 5, seed = 0)
  ro <- ComplexHeatmap::row_order(ht)

  cluster_sizes <- lengths(ro)
  expect_true(
    all(cluster_sizes > 0),
    info = paste("Empty clusters found. Sizes:", paste(cluster_sizes, collapse = ", "))
  )
})

test_that("B5: HeatmapList (draw output) has no @matrix slot", {
  skip_if_not_installed("ComplexHeatmap")
  # Documents the root cause of the original bug: kmeans_gene_lists used
  # shiny_env$ht@matrix, but draw() returns HeatmapList which has no @matrix.
  # The fix stores ht_matrix separately at render time.
  ht <- run_kmeans_heatmap(k = 4, seed = 0)

  expect_s4_class(ht, "HeatmapList")
  expect_false(
    "matrix" %in% slotNames(ht),
    info = "HeatmapList must NOT have a @matrix slot — confirms the fix is needed"
  )
})


# ===========================================================================
# SECTION C — set.seed() reproducibility
#
# The fix moves set.seed() to immediately before ComplexHeatmap::draw().
# Same seed must produce identical cluster assignments; different seeds
# are not guaranteed to differ but usually will.
# ===========================================================================

test_that("C1: same seed produces identical row_order across two runs", {
  ht_a <- run_kmeans_heatmap(k = 5, seed = 99)
  ht_b <- run_kmeans_heatmap(k = 5, seed = 99)

  ro_a <- ComplexHeatmap::row_order(ht_a)
  ro_b <- ComplexHeatmap::row_order(ht_b)

  # Sort by cluster name for a stable comparison
  expect_identical(ro_a[order(names(ro_a))], ro_b[order(names(ro_b))],
    info = "Same seed must produce identical cluster assignments")
})

test_that("C2: different seeds typically produce different row_order (stochastic)", {
  # This test documents expected behavior rather than asserting a hard guarantee.
  # If it fails, k-means may have converged to the same optimum — that's fine.
  ht_0 <- run_kmeans_heatmap(k = 8, seed = 0)
  ht_1 <- run_kmeans_heatmap(k = 8, seed = 1)

  ro_0 <- ComplexHeatmap::row_order(ht_0)
  ro_1 <- ComplexHeatmap::row_order(ht_1)

  sizes_0 <- sort(lengths(ro_0))
  sizes_1 <- sort(lengths(ro_1))

  # Sizes will differ when seeds produce different partitions.
  # If they happen to be identical, skip rather than fail.
  if (identical(sizes_0, sizes_1)) {
    skip("Seeds 0 and 1 converged to the same partition — stochastic skip")
  }
  expect_false(identical(sizes_0, sizes_1))
})


# ===========================================================================
# SECTION D — n_genes guard (transient-state crash prevention)
#
# heatmap_main_object() now requires nrow(heatmap_data()) > k_clusters.
# The guard lives in the module reactive (not testable directly), but we
# can confirm process_heatmap_data() never returns fewer than 10 rows.
# ===========================================================================

test_that("D1: process_heatmap_data never returns 0 rows (n_genes clamped to 10)", {
  result <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 0,        # invalid — should clamp to 10
    gene_centering = TRUE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  expect_gte(nrow(result), 10)
})

test_that("D2: heatmap_main() runs without error when k < nrow(data)", {
  small_data <- process_heatmap_data(
    data           = test_mat,
    n_genes_max    = 30,
    gene_centering = TRUE,
    gene_normalize = FALSE,
    sample_centering = FALSE,
    sample_normalize = FALSE,
    all_gene_names = dummy_gene_names,
    select_gene_id = "User_ID"
  )
  # k=5 < nrow=30, so this must succeed
  expect_no_error(run_kmeans_heatmap(k = 5, seed = 0, data = small_data))
})
