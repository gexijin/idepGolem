# Tests for normalized expression data (format = 2)
#
# Fixtures:
#   testdata/normalized_expr_testing/normalized_expr.csv  — 200 genes x 9 samples, log2(counts + 1),
#                                   columns: A_1, A_2, A_3, B_1, B_2, B_3, C_1, C_2, C_3
#   testdata/normalized_expr_testing/sample_info_normalized.csv — samples-as-rows, one column: condition (A/B/C)


# For manual testing, see manual_testing_checklist.md in the testdata/normalized_expr_testing directory.
# These tests are not intended to be exhaustive, but they guard the critical code paths for normalized expression data.


# ---------------------------------------------------------------------------
# Load fixtures once
# ---------------------------------------------------------------------------

norm_path  <- testthat::test_path("testdata", "normalized_expr_testing", "normalized_expr.csv")
sinfo_path <- testthat::test_path("testdata", "normalized_expr_testing", "sample_info_normalized.csv")

norm_expr  <- as.matrix(read.csv(norm_path,  row.names = 1, check.names = FALSE))
sample_info <- read.csv(sinfo_path, row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)

all_sample_names <- colnames(norm_expr)   # A_1 ... C_3

# ===========================================================================
# SECTION A — find_contrast_samples() with format=2 + design file
#
# Guards the bug fixed in fct_analysis_random.R: normalized data (format=2)
# must use the limma/else branch to resolve contrast samples from sample_info,
# not fall through to the default all-samples path.
# ===========================================================================

# A1. Simple contrast "A-B", single factor, format=2
test_that("A1: find_contrast_samples returns 6 indices for A-B (format=2)", {
  result <- find_contrast_samples(
    select_contrast        = "A-B",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: A vs. B",
    counts_deg_method      = 1,   # always provide; NULL causes if-guard to error
    data_file_format       = 2
  )
  expect_length(result, 6)
  expect_true(all(result %in% 1:6))
})

# A2. Reversed contrast "B-A" — same 6 indices as A1 (set equality).
#     The model_comprions string is stated as "condition: B vs. A" so that
#     the parsed comparison "B-A" matches select_contrast "B-A".
test_that("A2: find_contrast_samples returns same 6 indices for reversed B-A (format=2)", {
  result_ab <- find_contrast_samples(
    select_contrast        = "A-B",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: A vs. B",
    counts_deg_method      = 1,
    data_file_format       = 2
  )
  result_ba <- find_contrast_samples(
    select_contrast        = "B-A",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: B vs. A",
    counts_deg_method      = 1,
    data_file_format       = 2
  )
  expect_setequal(result_ab, result_ba)
  expect_length(result_ba, 6)
  expect_true(all(result_ba %in% 1:6))
})

# A3. Contrast string "A-B" not present in select_model_comprions ("B-C") —
#     falls back to all 9 indices.
test_that("A3: find_contrast_samples falls back to all samples when contrast not in comparisons", {
  result <- find_contrast_samples(
    select_contrast        = "A-B",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: B vs. C",
    counts_deg_method      = 1,
    data_file_format       = 2
  )
  expect_length(result, 9)
  expect_equal(sort(result), seq_along(all_sample_names))
})

# A4. sample_info = NULL (no design file), 2-group column names A_1..B_3.
#     detect_groups() strips trailing numbers/underscores: A_1 -> A, B_1 -> B.
#     Contrast "A-B" matches groups A and B → indices 1:6.
test_that("A4: find_contrast_samples uses column-name groups when sample_info is NULL", {
  result <- find_contrast_samples(
    select_contrast  = "A-B",
    all_sample_names = all_sample_names,
    sample_info      = NULL,
    data_file_format = 2
  )
  expect_length(result, 6)
  expect_true(all(result %in% 1:6))
})

# A5. format=1 + limma-trend (counts_deg_method = 1), same design as A1.
#     The fixed else-branch (not DESeq2) must still resolve the correct 6 indices.
test_that("A5: find_contrast_samples returns 6 indices for format=1 limma-trend", {
  result <- find_contrast_samples(
    select_contrast        = "A-B",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: A vs. B",
    counts_deg_method      = 1,
    data_file_format       = 1
  )
  expect_length(result, 6)
  expect_true(all(result %in% 1:6))
})

# A6. Interaction term contrast "I:condition_batch" — always uses all samples.
test_that("A6: find_contrast_samples returns all 9 indices for interaction term contrast", {
  result <- find_contrast_samples(
    select_contrast        = "I:condition_batch",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: A vs. B",
    counts_deg_method      = 1,
    data_file_format       = 2
  )
  expect_length(result, 9)
  expect_equal(sort(result), seq_along(all_sample_names))
})

# A7. format=1 + DESeq2 (counts_deg_method = 3) — regression guard that the
#     DESeq2 code path was not broken by the format=2 fix.
test_that("A7: find_contrast_samples returns 6 indices for format=1 DESeq2", {
  result <- find_contrast_samples(
    select_contrast        = "A-B",
    all_sample_names       = all_sample_names,
    sample_info            = sample_info,
    select_factors_model   = "condition",
    select_model_comprions = "condition: A vs. B",
    counts_deg_method      = 3,
    data_file_format       = 1
  )
  expect_length(result, 6)
  expect_true(all(result %in% 1:6))
})

# ===========================================================================
# SECTION B — deg_limma() return structure for format=2
#
# Confirms the list shape, column names, and error-path behaviour of
# deg_limma when called with normalised expression data.
# ===========================================================================

# Subset fixtures for B tests
norm_expr_ab  <- norm_expr[, 1:6]          # A and B columns only (6 samples)
norm_expr_1pg <- norm_expr[, c(1, 4, 7)]  # one sample per group A/B/C

# Block-annotated sample_info for B3:
# "rep" is balanced across conditions, so duplicateCorrelation should succeed.
sinfo_block <- data.frame(
  condition = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
  rep       = c("rep1", "rep2", "rep3", "rep1", "rep2", "rep3",
                "rep1", "rep2", "rep3"),
  row.names = colnames(norm_expr),
  stringsAsFactors = FALSE
)

# B1. Single comparison A vs. B, format=2 (2-group path, no design file)
# -------------------------------------------------------------------------
# With only 2 unique groups and model_factors=NULL, deg_limma takes the
# 2-group limma path, which explicitly renames "logFC" -> "log2FC".
test_that("B1: deg_limma returns correct list structure for format=2, A vs. B", {
  result <- deg_limma(
    processed_data       = norm_expr_ab,
    max_p_limma          = 0.99,
    min_fc_limma         = 0,
    raw_counts           = norm_expr_ab,
    counts_deg_method    = 1,
    prior_counts         = 1,
    data_file_format     = 2,
    selected_comparisons = NULL,
    sample_info          = NULL,
    model_factors        = NULL,
    block_factor         = NULL,
    reference_levels     = NULL
  )
  expect_true(is.list(result))
  expect_true(all(c("results", "comparisons", "top_genes", "expr") %in% names(result)))
  expect_equal(ncol(result$top_genes[[1]]), 2)
  expect_equal(colnames(result$top_genes[[1]]), c("log2FC", "adj.P.Val"))
  expect_true(nrow(result$results) > 0)
})

# B2. Two comparisons selected (A vs. B and B vs. C)
# -------------------------------------------------------------------------
# Uses model_factors + design file -> 3+ group path -> two named top_genes.
test_that("B2: deg_limma returns two named comparisons for format=2 multi-comparison", {
  result <- deg_limma(
    processed_data       = norm_expr,
    max_p_limma          = 0.99,
    min_fc_limma         = 0,
    raw_counts           = norm_expr,
    counts_deg_method    = 1,
    prior_counts         = 1,
    data_file_format     = 2,
    selected_comparisons = c("condition: A vs. B", "condition: B vs. C"),
    sample_info          = sample_info,
    model_factors        = c("condition"),
    block_factor         = NULL,
    reference_levels     = NULL
  )
  expect_true(is.list(result))
  expect_equal(length(result$comparisons), 2)
  expect_equal(length(result$top_genes), 2)
  expect_equal(names(result$top_genes), result$comparisons)
})

# B3. format=2 with a balanced block factor (paired/replicate design)
# -------------------------------------------------------------------------
# "rep" factor is balanced: rep1/rep2/rep3 appear once in each condition,
# so duplicateCorrelation should not warn about confounding.
test_that("B3: deg_limma with block factor returns valid list for format=2", {
  result <- deg_limma(
    processed_data       = norm_expr,
    max_p_limma          = 0.99,
    min_fc_limma         = 0,
    raw_counts           = norm_expr,
    counts_deg_method    = 1,
    prior_counts         = 1,
    data_file_format     = 2,
    selected_comparisons = c("condition: A vs. B"),
    sample_info          = sinfo_block,
    model_factors        = c("condition"),
    block_factor         = "rep",
    reference_levels     = NULL
  )
  expect_true(is.list(result))
  expect_false(is.character(result))
})

# B4. format=2 with only 1 replicate per group
# -------------------------------------------------------------------------
# reps check: sum(reps[,1] >= 2) == 0 < 2 -> returns named list with
# exp_type error, never throws.
test_that("B4: deg_limma returns exp_type error (not a throw) for 1-rep-per-group format=2", {
  result <- deg_limma(
    processed_data       = norm_expr_1pg,
    max_p_limma          = 0.99,
    min_fc_limma         = 0,
    raw_counts           = norm_expr_1pg,
    counts_deg_method    = 1,
    prior_counts         = 1,
    data_file_format     = 2,
    selected_comparisons = NULL,
    sample_info          = NULL,
    model_factors        = NULL,
    block_factor         = NULL,
    reference_levels     = NULL
  )
  expect_true(is.list(result))
  expect_false(is.null(result$exp_type))
  expect_true(is.character(result$exp_type))
})

# ===========================================================================
# SECTION C — deg_information() for format=2
#
# deg_information joins limma results with gene metadata. Tests use a
# small synthetic limma_value so they run without the full app stack.
#
# Return structure: list(degs_data, limma_value$Results)
#   result[[1]] = data frame of gene-level stats
#   result[[2]] = decideTests matrix (may be NULL for mocks)
# ===========================================================================

# Build a small mock limma_value (50 genes, 1 comparison "A-B")
# Mimics deg_limma output for the 2-group path (format=2, no baseMean).
local({
  n_mock <- 50L
  gids   <- paste0("gene", sprintf("%04d", seq_len(n_mock)))

  top_ab <- data.frame(
    log2FC    = rnorm(n_mock, mean = 1),
    adj.P.Val = runif(n_mock),
    row.names = gids,
    stringsAsFactors = FALSE
  )
  top_bc <- data.frame(
    log2FC    = rnorm(n_mock, mean = -0.5),
    adj.P.Val = runif(n_mock),
    row.names = gids,
    stringsAsFactors = FALSE
  )

  limma_mock_1 <<- list(
    results     = matrix(
      sample(c(-1L, 0L, 1L), n_mock, replace = TRUE),
      nrow = n_mock, ncol = 1,
      dimnames = list(gids, "A-B")
    ),
    comparisons = "A-B",
    top_genes   = setNames(list(top_ab), "A-B"),
    baseMean    = NULL,
    exp_type    = "2 sample groups."
  )

  limma_mock_2 <<- list(
    results     = matrix(
      sample(c(-1L, 0L, 1L), n_mock * 2, replace = TRUE),
      nrow = n_mock, ncol = 2,
      dimnames = list(gids, c("A-B", "B-C"))
    ),
    comparisons = c("A-B", "B-C"),
    top_genes   = setNames(list(top_ab, top_bc), c("A-B", "B-C")),
    baseMean    = NULL,
    exp_type    = "3 sample groups detected."
  )

  # gene_names: include User_ID, ensembl_ID, symbol so that
  # deg_information's dplyr::relocate() calls succeed for both paths.
  gene_names_mock <<- data.frame(
    User_ID    = gids,
    ensembl_ID = gids,
    symbol     = paste0("SYM", seq_len(n_mock)),
    stringsAsFactors = FALSE
  )

  # Minimal processed_data matrix matching the mock gene IDs
  set.seed(42)
  processed_mock <<- matrix(
    rnorm(n_mock * 6),
    nrow = n_mock, ncol = 6,
    dimnames = list(
      gids,
      c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    )
  )
})

# C1. no_id_conversion = TRUE (gene symbols as row names, no Ensembl lookup)
test_that("C1: deg_information with no_id_conversion=TRUE returns User_ID and LFC/adjPval columns", {
  result <- deg_information(
    limma_value      = limma_mock_1,
    gene_names       = gene_names_mock,
    processed_data   = processed_mock,
    no_id_conversion = TRUE
  )
  degs <- result[[1]]
  expect_true(is.data.frame(degs))
  expect_true(nrow(degs) > 0)
  expect_true("User_ID" %in% colnames(degs))
  # LFC and adjPval columns are named "<comparison>_log2FC" / "<comparison>_adjPval"
  expect_true("A-B_log2FC"  %in% colnames(degs))
  expect_true("A-B_adjPval" %in% colnames(degs))
})

# C2. no_id_conversion = FALSE (Ensembl ID path)
test_that("C2: deg_information with no_id_conversion=FALSE returns ensembl_ID and symbol columns", {
  result <- deg_information(
    limma_value      = limma_mock_1,
    gene_names       = gene_names_mock,
    processed_data   = processed_mock,
    no_id_conversion = FALSE
  )
  degs <- result[[1]]
  expect_true(is.data.frame(degs))
  expect_true(nrow(degs) > 0)
  expect_true("ensembl_ID" %in% colnames(degs))
  expect_true("symbol" %in% colnames(degs))
})

# C3. Multi-comparison result — both LFC/adjPval pairs present
test_that("C3: deg_information with 2 comparisons includes both LFC and adjPval column pairs", {
  result <- deg_information(
    limma_value      = limma_mock_2,
    gene_names       = gene_names_mock,
    processed_data   = processed_mock,
    no_id_conversion = TRUE
  )
  degs <- result[[1]]
  expect_true("A-B_log2FC"  %in% colnames(degs))
  expect_true("A-B_adjPval" %in% colnames(degs))
  expect_true("B-C_log2FC"  %in% colnames(degs))
  expect_true("B-C_adjPval" %in% colnames(degs))
})

# C4. limma result (baseMean = NULL) — no baseMean column in output
test_that("C4: deg_information omits baseMean column when limma_value$baseMean is NULL", {
  # limma_mock_1 already has baseMean = NULL (limma, not DESeq2)
  result <- deg_information(
    limma_value      = limma_mock_1,
    gene_names       = gene_names_mock,
    processed_data   = processed_mock,
    no_id_conversion = TRUE
  )
  degs <- result[[1]]
  expect_false("baseMean" %in% colnames(degs))
})

# ===========================================================================
# Shared fixtures for Sections D, E, F
#
# A self-contained 50-gene × 9-sample dataset so these tests do not depend
# on row-name alignment between norm_expr and the Section C mocks.
# ===========================================================================

local({
  set.seed(42)
  n_de  <- 50L
  genes <- paste0("gene", sprintf("%04d", seq_len(n_de)))

  de_mat <<- matrix(
    rnorm(n_de * 9),
    nrow     = n_de,
    ncol     = 9,
    dimnames = list(
      genes,
      c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3", "C_1", "C_2", "C_3")
    )
  )

  # top_genes entry: log2FC and adj.P.Val for every gene
  de_top <<- data.frame(
    log2FC    = rnorm(n_de),
    adj.P.Val = runif(n_de, 0.001, 0.5),
    row.names = genes,
    stringsAsFactors = FALSE
  )
  de_top_list <<- setNames(list(de_top), "A-B")

  # results matrix: force first 10 up, next 10 down, rest zero
  de_results <<- matrix(
    0L,
    nrow     = n_de,
    ncol     = 1,
    dimnames = list(genes, "A-B")
  )
  de_results[1:10,  1] <<-  1L
  de_results[11:20, 1] <<- -1L

  de_limma <<- list(
    results     = de_results,
    comparisons = "A-B",
    top_genes   = de_top_list,
    baseMean    = NULL,
    exp_type    = "2 sample groups."
  )
})

# ===========================================================================
# SECTION D — volcano_data() sample subsetting
#
# Guards that Average is computed over contrast_samples only,
# not all columns of processed_data.
# ===========================================================================

# D1. contrast_samples = 1:6 → Average equals rowMeans(mat[, 1:6])
test_that("D1: volcano_data Average equals rowMeans of contrast columns only", {
  result <- volcano_data(
    select_contrast  = "A-B",
    comparisons      = "A-B",
    top_genes        = de_top_list,
    limma_p_val      = 0.99,
    limma_fc         = 1,
    processed_data   = de_mat,
    contrast_samples = 1:6,
    all_gene_names   = NULL,
    select_gene_id   = NULL
  )
  expect_false(is.null(result))

  ref6 <- rowMeans(de_mat[, 1:6])
  ref6_df <- data.frame(Row.names = rownames(de_mat), ref6 = ref6,
                        stringsAsFactors = FALSE)
  cmp <- merge(result$data, ref6_df, by = "Row.names")
  expect_equal(cmp$Average, cmp$ref6, tolerance = 1e-10)

  # Must differ from full-matrix average
  ref9 <- rowMeans(de_mat[, 1:9])
  ref9_df <- data.frame(Row.names = rownames(de_mat), ref9 = ref9,
                        stringsAsFactors = FALSE)
  cmp9 <- merge(result$data, ref9_df, by = "Row.names")
  expect_false(isTRUE(all.equal(cmp9$Average, cmp9$ref9)))
})

# D2. contrast_samples = 1:6 vs. 1:9 → Average values differ between calls
test_that("D2: volcano_data Average differs when contrast_samples differ", {
  res6 <- volcano_data(
    select_contrast  = "A-B",
    comparisons      = "A-B",
    top_genes        = de_top_list,
    limma_p_val      = 0.99,
    limma_fc         = 1,
    processed_data   = de_mat,
    contrast_samples = 1:6,
    all_gene_names   = NULL,
    select_gene_id   = NULL
  )
  res9 <- volcano_data(
    select_contrast  = "A-B",
    comparisons      = "A-B",
    top_genes        = de_top_list,
    limma_p_val      = 0.99,
    limma_fc         = 1,
    processed_data   = de_mat,
    contrast_samples = 1:9,
    all_gene_names   = NULL,
    select_gene_id   = NULL
  )
  # Averages must differ between the two calls
  avg6 <- res6$data$Average[order(res6$data$Row.names)]
  avg9 <- res9$data$Average[order(res9$data$Row.names)]
  expect_false(isTRUE(all.equal(avg6, avg9)))
})

# D3. Up/Down labels respect p-value and fold-change thresholds
test_that("D3: volcano_data Up/Down labeling respects limma_p_val and limma_fc thresholds", {
  p_thr <- 0.3
  fc_thr <- 1.5
  result <- volcano_data(
    select_contrast  = "A-B",
    comparisons      = "A-B",
    top_genes        = de_top_list,
    limma_p_val      = p_thr,
    limma_fc         = fc_thr,
    processed_data   = de_mat,
    contrast_samples = 1:6,
    all_gene_names   = NULL,
    select_gene_id   = NULL
  )
  d <- result$data
  up_rows   <- d[d$upOrDown == "Up",   ]
  down_rows <- d[d$upOrDown == "Down", ]
  if (nrow(up_rows) > 0) {
    expect_true(all(up_rows$FDR  <= p_thr))
    expect_true(all(up_rows$Fold >= log2(fc_thr)))
  }
  if (nrow(down_rows) > 0) {
    expect_true(all(down_rows$FDR  <= p_thr))
    expect_true(all(down_rows$Fold <= -log2(fc_thr)))
  }
})

# ===========================================================================
# SECTION E — deg_heat_data() sample subsetting
#
# Guards that the heatmap matrix uses only contrast_samples columns.
# ===========================================================================

# E1. Returns a matrix with exactly contrast_samples columns and a valid bar
test_that("E1: deg_heat_data returns only contrast-sample columns", {
  result <- deg_heat_data(
    limma            = de_limma,
    select_contrast  = "A-B",
    processed_data   = de_mat,
    contrast_samples = 1:6
  )
  expect_false(is.null(result))
  expect_equal(ncol(result$genes), 6)
  expect_false(is.null(result$bar))
  expect_true(all(result$bar %in% c(-1, 1)))
})

# E2. No significant genes → returns NULL without error
test_that("E2: deg_heat_data returns NULL when no significant genes", {
  de_limma_zero <- de_limma
  de_limma_zero$results <- matrix(
    0L,
    nrow     = nrow(de_results),
    ncol     = 1,
    dimnames = list(rownames(de_results), "A-B")
  )
  expect_null(
    deg_heat_data(
      limma            = de_limma_zero,
      select_contrast  = "A-B",
      processed_data   = de_mat,
      contrast_samples = 1:6
    )
  )
})

# ===========================================================================
# SECTION F — build_deg_gene_table() (gene list download content)
#
# build_deg_gene_table() is the pure-function extraction of the
# deg_gene_table_data reactive; it can be tested without Shiny.
# ===========================================================================

# F1. Meta and display structure for format=2 (no gene_names / gene_info)
test_that("F1: build_deg_gene_table returns correct meta and display columns", {
  result <- build_deg_gene_table(
    top_list              = de_top_list,
    select_contrast       = "A-B",
    limma_p_val           = 0.99,
    limma_fc              = 1,
    gene_direction_filter = "Up-regulated genes",
    gene_names            = NULL,
    gene_info             = NULL
  )

  expect_true(nrow(result$meta) > 0)

  # meta column names
  expect_true(all(
    c("ensembl_ID", "symbol", "entrezgene_id",
      "log2FC", "Adjusted_P_value", "description") %in% colnames(result$meta)
  ))

  # display column names
  expect_true(all(
    c("ID", "Ensembl ID", "log2 FC", "Adj. Pval", "Description") %in%
      colnames(result$display)
  ))

  # order_dir is one of the two valid values
  expect_true(result$order_dir %in% c("asc", "desc"))
})

# F2. meta is a valid writable/readable CSV
test_that("F2: build_deg_gene_table meta round-trips through write.csv / read.csv", {
  result <- build_deg_gene_table(
    top_list              = de_top_list,
    select_contrast       = "A-B",
    limma_p_val           = 0.99,
    limma_fc              = 1,
    gene_direction_filter = "Up-regulated genes",
    gene_names            = NULL,
    gene_info             = NULL
  )

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error(write.csv(result$meta, tmp, row.names = FALSE))

  back <- read.csv(tmp, stringsAsFactors = FALSE)
  expect_equal(
    sort(colnames(back)),
    sort(colnames(result$meta))
  )
})
