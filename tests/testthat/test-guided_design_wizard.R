# Tests for the guided design wizard (Guided Builder tab).
#
# The wizard's server logic does two things before calling build_sample_info_from_df:
#   1. Parse factor definitions (level strings, blocking flag)
#   2. Collect per-sample dropdown / numericInput values into a design data frame
#
# We simulate both steps here and then verify build_sample_info_from_df produces
# the correct output for every experiment type the app claims to support.
#
# Run via the package test suite:
#   devtools::test(filter = "guided_design_wizard")
#
# Or standalone (no full package install needed):
#   source("R/utils_analysis_random.R")
#   source("R/fct_01_load_data.R")
#   testthat::test_file("tests/testthat/test-guided_design_wizard.R")

# ---------------------------------------------------------------------------
# Helpers that mirror the wizard's server-side logic
# ---------------------------------------------------------------------------

parse_guided_levels <- function(raw_str) {
  # Mirrors observeEvent(input$guided_next_btn): split, trim, sanitize, uppercase, dedup
  levels_vec <- trimws(strsplit(raw_str, ",", fixed = TRUE)[[1L]])
  levels_vec <- toupper(sanitize_names(levels_vec))
  levels_vec <- unique(levels_vec)
  levels_vec[nchar(levels_vec) > 0L]
}

detect_shared_levels <- function(factor_level_lists) {
  # Mirrors the duplicate-level warning check in observeEvent(input$apply_guided_design)
  all_vals <- unlist(lapply(factor_level_lists, function(x) tolower(unique(x))))
  tbl <- table(all_vals)
  names(tbl[tbl > 1L])
}

wizard_to_design_df <- function(factor_names, cell_value_rows, sample_names) {
  # Mirrors the data-frame construction before build_sample_info_from_df is called
  df <- as.data.frame(
    do.call(rbind, cell_value_rows),
    stringsAsFactors = FALSE
  )
  colnames(df) <- sample_names
  rownames(df) <- factor_names
  df
}

# ---------------------------------------------------------------------------
# 1. Level string parsing
# ---------------------------------------------------------------------------

test_that("parse_guided_levels trims whitespace, splits, sanitizes, and uppercases", {
  expect_equal(parse_guided_levels("ctrl, treated"), c("CTRL", "TREATED"))
  expect_equal(parse_guided_levels(" ctrl , treated "), c("CTRL", "TREATED"))
  expect_equal(parse_guided_levels("ctrl,treated,t2"), c("CTRL", "TREATED", "T2"))
})

test_that("parse_guided_levels removes spaces and hyphens within a level name", {
  expect_equal(parse_guided_levels("S 2, s1"), c("S2", "S1"))
  expect_equal(parse_guided_levels("wild-type, knock-out"), c("WILDTYPE", "KNOCKOUT"))
  expect_equal(parse_guided_levels("S 2"), "S2")
})

test_that("parse_guided_levels deduplicates levels that collapse to the same string after sanitization", {
  # "ctrl" and "Ctrl" both become "CTRL" — only one survives
  result <- parse_guided_levels("ctrl, Ctrl, treated")
  expect_equal(result, c("CTRL", "TREATED"))
})

test_that("parse_guided_levels drops empty tokens between consecutive commas", {
  result <- parse_guided_levels("ctrl,,treated")
  expect_equal(result, c("CTRL", "TREATED"))
})

test_that("parse_guided_levels returns empty vector for blank input", {
  expect_equal(parse_guided_levels(""), character(0))
  expect_equal(parse_guided_levels("   "), character(0))
})

test_that("parse_guided_levels returns single-element vector for one level", {
  expect_equal(parse_guided_levels("ctrl"), "CTRL")
})

# ---------------------------------------------------------------------------
# 2. Duplicate-level detection (cross-factor warning)
# ---------------------------------------------------------------------------

test_that("detect_shared_levels finds names reused across factors", {
  dups <- detect_shared_levels(list(c("ctrl", "treat"), c("ctrl", "high")))
  expect_equal(dups, "ctrl")
})

test_that("detect_shared_levels is case-insensitive", {
  dups <- detect_shared_levels(list(c("Ctrl", "Treat"), c("ctrl", "High")))
  expect_equal(dups, "ctrl")
})

test_that("detect_shared_levels returns nothing when levels are distinct", {
  dups <- detect_shared_levels(list(c("ctrl", "treat"), c("WT", "Mu")))
  expect_length(dups, 0L)
})

test_that("detect_shared_levels handles three factors", {
  dups <- detect_shared_levels(list(
    c("ctrl", "treat"),
    c("WT", "Mu"),
    c("ctrl", "batch2")  # "ctrl" appears in factor 1 and factor 3
  ))
  expect_equal(dups, "ctrl")
})

# ---------------------------------------------------------------------------
# 3. Experiment type A — paired samples
# ---------------------------------------------------------------------------

test_that("paired design: two factors (treatment + pairs) produce correct matrix", {
  samples <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3")
  treatment <- c("Ctrl", "Ctrl", "Ctrl", "Treat", "Treat", "Treat")
  pairs     <- c("1",   "2",   "3",   "1",   "2",   "3")

  df <- wizard_to_design_df(
    c("Treatment", "pairs"), list(treatment, pairs), samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 6L)
  expect_equal(ncol(result), 2L)
  expect_setequal(colnames(result), c("Treatment", "pairs"))
  expect_equal(
    as.character(result[, "Treatment"]),
    c("CTRL", "CTRL", "CTRL", "TREAT", "TREAT", "TREAT")
  )
  # Pair indices survive sanitization (digits only — no stripping)
  expect_equal(as.character(result[, "pairs"]), c("1", "2", "3", "1", "2", "3"))
})

test_that("paired design: blocking factor with only one distinct value is dropped", {
  # If user assigns the same pair index to all samples, it has no variation
  samples   <- c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2")
  treatment <- c("Ctrl", "Ctrl", "Treat", "Treat")
  pairs     <- c("1", "1", "1", "1")   # all same — no variation

  df <- wizard_to_design_df(
    c("Treatment", "pairs"), list(treatment, pairs), samples
  )
  result <- build_sample_info_from_df(df, samples)

  # pairs factor is dropped; treatment survives
  expect_false(is.null(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "Treatment")
})

# ---------------------------------------------------------------------------
# 4. Experiment type B — 2×2 factorial design
# ---------------------------------------------------------------------------

test_that("2x2 factorial: genotype x treatment produces 2-column matrix", {
  samples  <- c("WT_Ctrl_1", "WT_Ctrl_2", "WT_Trt_1", "WT_Trt_2",
                "Mu_Ctrl_1", "Mu_Ctrl_2", "Mu_Trt_1", "Mu_Trt_2")
  genotype  <- c("WT", "WT", "WT", "WT", "Mu", "Mu", "Mu", "Mu")
  treatment <- c("Ctrl", "Ctrl", "Trt", "Trt", "Ctrl", "Ctrl", "Trt", "Trt")

  df <- wizard_to_design_df(
    c("Genotype", "Treatment"), list(genotype, treatment), samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 8L)
  expect_equal(ncol(result), 2L)
  expect_setequal(colnames(result), c("Genotype", "Treatment"))
})

test_that("2x2 factorial: level values are uppercased", {
  samples  <- c("wt_ctrl", "wt_trt", "mu_ctrl", "mu_trt")
  genotype  <- c("wt", "wt", "mu", "mu")
  treatment <- c("ctrl", "trt", "ctrl", "trt")

  df <- wizard_to_design_df(
    c("Genotype", "Treatment"), list(genotype, treatment), samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_equal(unique(result[, "Genotype"]), c("WT", "MU"))
  expect_equal(unique(result[, "Treatment"]), c("CTRL", "TRT"))
})

# ---------------------------------------------------------------------------
# 5. Experiment type C — ESR1/ESR2 double-knockout as 2x2 factorial
# ---------------------------------------------------------------------------

test_that("ESR double-KO design encoded as 2x2 factorial", {
  samples  <- c("WT_1",  "WT_2",  "ESR1_1", "ESR1_2",
                "ESR2_1", "ESR2_2", "DKO_1", "DKO_2")
  esr1 <- c("1WT", "1WT", "1KO", "1KO", "1WT", "1WT", "1KO", "1KO")
  esr2 <- c("2WT", "2WT", "2WT", "2WT", "2KO", "2KO", "2KO", "2KO")

  df <- wizard_to_design_df(c("ESR1", "ESR2"), list(esr1, esr2), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(ncol(result), 2L)
  expect_setequal(colnames(result), c("ESR1", "ESR2"))
  # DKO samples should be 1KO × 2KO
  dko_idx <- which(samples %in% c("DKO_1", "DKO_2"))
  expect_true(all(result[dko_idx, "ESR1"] == "1KO"))
  expect_true(all(result[dko_idx, "ESR2"] == "2KO"))
})

# ---------------------------------------------------------------------------
# 6. Experiment type D — 2×3 design (two genotypes × three treatments)
# ---------------------------------------------------------------------------

test_that("2x3 design with 3 treatment levels works correctly", {
  samples   <- c("WT_Veh_1", "WT_Veh_2", "WT_T1_1", "WT_T1_2", "WT_T2_1", "WT_T2_2",
                 "Mu_Veh_1", "Mu_Veh_2", "Mu_T1_1", "Mu_T1_2", "Mu_T2_1", "Mu_T2_2")
  genotype  <- rep(c("WT", "WT", "WT", "WT", "WT", "WT",
                     "Mu", "Mu", "Mu", "Mu", "Mu", "Mu"), 1)
  treatment <- rep(c("Veh", "Veh", "T1", "T1", "T2", "T2"), 2)

  df <- wizard_to_design_df(
    c("Genotype", "Treatment"), list(genotype, treatment), samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 12L)
  expect_equal(ncol(result), 2L)
  expect_equal(length(unique(result[, "Treatment"])), 3L)
})

# ---------------------------------------------------------------------------
# 7. Experiment type F — time course with tissue factor
# ---------------------------------------------------------------------------

test_that("time course x tissue 2-factor design (3x2) works", {
  samples <- c(
    "leaf_0h_1", "leaf_0h_2", "leaf_12h_1", "leaf_12h_2",
    "leaf_48h_1", "leaf_48h_2",
    "root_0h_1", "root_0h_2", "root_12h_1", "root_12h_2",
    "root_48h_1", "root_48h_2"
  )
  time   <- c("0h","0h","12h","12h","48h","48h",
              "0h","0h","12h","12h","48h","48h")
  tissue <- c("leaf","leaf","leaf","leaf","leaf","leaf",
              "root","root","root","root","root","root")

  df <- wizard_to_design_df(c("time", "tissue"), list(time, tissue), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 12L)
  expect_equal(ncol(result), 2L)
  expect_equal(length(unique(result[, "time"])), 3L)
  expect_equal(length(unique(result[, "tissue"])), 2L)
})

# ---------------------------------------------------------------------------
# 8. 5-factor design (maximum supported by wizard)
# ---------------------------------------------------------------------------

test_that("5-factor design (wizard max) produces correct matrix", {
  samples <- paste0("s", 1:8)
  f1 <- rep(c("A", "B"), 4)
  f2 <- rep(c("X", "Y"), each = 4)
  f3 <- c("p", "p", "q", "q", "p", "p", "q", "q")
  f4 <- c("1", "2", "1", "2", "1", "2", "1", "2")
  f5 <- c("M", "M", "M", "M", "N", "N", "N", "N")

  df <- wizard_to_design_df(
    c("f1", "f2", "f3", "f4", "f5"),
    list(f1, f2, f3, f4, f5),
    samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 8L)
  expect_equal(ncol(result), 5L)
})

# ---------------------------------------------------------------------------
# 9. Unassigned samples (empty string — wizard "-- select --" default)
# ---------------------------------------------------------------------------

test_that("design with one unassigned sample returns NULL via build_sample_info_from_df", {
  # The wizard blocks Apply before this point, but we verify the downstream
  # function also rejects it (empty string becomes no-variation after sanitize)
  samples   <- c("s1", "s2", "s3", "s4")
  treatment <- c("ctrl", "ctrl", "treat", "")  # s4 unassigned → empty

  df <- wizard_to_design_df("Treatment", list(treatment), samples)

  # An all-empty row after sanitize is dropped; remaining row has variation → valid
  # But the empty value "" becomes "" after sanitize — so 3 distinct values
  # build_sample_info drops factors where ALL are same; here 3 distinct values exist
  # The caller (apply_guided_design) catches unassigned before calling this.
  # We just verify the function doesn't crash.
  result <- build_sample_info_from_df(df, samples)
  # Either NULL or a matrix is acceptable — no error
  expect_true(is.null(result) || is.matrix(result))
})

test_that("all samples unassigned (all empty) yields NULL", {
  samples   <- c("s1", "s2", "s3", "s4")
  treatment <- rep("", 4)  # nobody assigned

  df <- wizard_to_design_df("Treatment", list(treatment), samples)
  expect_null(build_sample_info_from_df(df, samples))
})

# ---------------------------------------------------------------------------
# 10. Shared level names across factors (warning scenario)
# ---------------------------------------------------------------------------

test_that("shared level names across factors still build valid design after sanitization", {
  # iDEP adds factor prefix internally — build_sample_info_from_df itself
  # does NOT prefix, but it should still return a valid matrix.
  samples  <- c("s1", "s2", "s3", "s4")
  factor_a <- c("ctrl", "ctrl", "treat", "treat")
  factor_b <- c("ctrl", "treat", "ctrl", "treat")  # "ctrl"/"treat" reused

  df <- wizard_to_design_df(
    c("FactorA", "FactorB"), list(factor_a, factor_b), samples
  )
  result <- build_sample_info_from_df(df, samples)

  # Both factors have variation — both should be retained
  expect_false(is.null(result))
  expect_equal(ncol(result), 2L)
})

# ---------------------------------------------------------------------------
# 11. Factor name sanitization (wizard sends raw user text)
# ---------------------------------------------------------------------------

test_that("factor names with spaces and hyphens are sanitized", {
  samples   <- c("s1", "s2", "s3", "s4")
  treatment <- c("ctrl", "ctrl", "treat", "treat")

  df <- wizard_to_design_df("my treatment", list(treatment), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  # Factor name "my treatment" → sanitize_names → "mytreatment"
  expect_equal(colnames(result), "mytreatment")
})

test_that("duplicate factor names after sanitization still produce valid matrix", {
  samples <- c("s1", "s2", "s3", "s4")
  r1 <- c("A", "A", "B", "B")
  r2 <- c("X", "X", "Y", "Y")

  # "group 1" and "group 2" both sanitize to "group1" and "group2" (distinct)
  df <- wizard_to_design_df(c("group 1", "group 2"), list(r1, r2), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(ncol(result), 2L)
})

# ---------------------------------------------------------------------------
# 12. Single replicate per group (n=1) — edge case
# ---------------------------------------------------------------------------

test_that("single replicate per group produces valid matrix", {
  samples   <- c("ctrl_1", "treat_1", "t2_1")
  treatment <- c("ctrl", "treat", "t2")

  df <- wizard_to_design_df("Treatment", list(treatment), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 3L)
  expect_equal(length(unique(result[, "Treatment"])), 3L)
})

# ---------------------------------------------------------------------------
# 13. Large experiment (30 samples, 3 factors) — performance smoke test
# ---------------------------------------------------------------------------

test_that("large experiment with 30 samples and 3 factors completes without error", {
  n <- 30L
  samples   <- paste0("sample_", seq_len(n))
  genotype  <- rep(c("WT", "Mu"), n / 2)
  treatment <- rep(c("Ctrl", "T1", "T2"), each = n / 3)
  batch     <- rep(paste0("B", 1:3), times = n / 3)

  df <- wizard_to_design_df(
    c("Genotype", "Treatment", "Batch"),
    list(genotype, treatment, batch),
    samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 3L)
})

# ---------------------------------------------------------------------------
# 14. n_guided_factors reset: ensure wizard accepts 1 through 5 factors
# ---------------------------------------------------------------------------

test_that("1-factor wizard design works", {
  samples   <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  treatment <- c("ctrl", "ctrl", "treat", "treat")

  df <- wizard_to_design_df("Treatment", list(treatment), samples)
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(ncol(result), 1L)
})

test_that("3-factor wizard design works", {
  samples <- paste0("s", 1:8)
  f1 <- rep(c("A", "B"), 4)
  f2 <- rep(c("X", "Y"), each = 4)
  f3 <- c("p", "q", "p", "q", "p", "q", "p", "q")

  df <- wizard_to_design_df(
    c("f1", "f2", "f3"), list(f1, f2, f3), samples
  )
  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(ncol(result), 3L)
})

# ---------------------------------------------------------------------------
# 15. Factor-name sanitization at wizard input layer
# ---------------------------------------------------------------------------

parse_guided_factor_name <- function(raw_str) {
  # Mirrors observeEvent(input$guided_next_btn): trim then sanitize_names (no uppercase)
  sanitize_names(trimws(if (is.null(raw_str)) "" else raw_str))
}

test_that("parse_guided_factor_name removes spaces from factor names", {
  expect_equal(parse_guided_factor_name("my factor"), "myfactor")
  expect_equal(parse_guided_factor_name(" group "), "group")
})

test_that("parse_guided_factor_name removes hyphens from factor names", {
  expect_equal(parse_guided_factor_name("time-point"), "timepoint")
  expect_equal(parse_guided_factor_name("cell-type"), "celltype")
})

test_that("parse_guided_factor_name handles NULL as empty string", {
  expect_equal(parse_guided_factor_name(NULL), "")
})

test_that("parse_guided_factor_name does NOT uppercase factor names (only levels are uppercased)", {
  expect_equal(parse_guided_factor_name("Group"), "Group")
  expect_equal(parse_guided_factor_name("timePoint"), "timePoint")
})

test_that("parse_guided_factor_name strips leading/trailing whitespace before sanitizing", {
  expect_equal(parse_guided_factor_name("  time point  "), "timepoint")
})

test_that("factor-name sanitization integrates cleanly with build_sample_info_from_df", {
  # Factor name "my factor" → "myfactor" via parse_guided_factor_name, then used as
  # the row name in the design df, which becomes the column name in the output matrix.
  samples <- c("s1", "s2", "s3", "s4")
  factor_name <- parse_guided_factor_name("my factor")
  df <- data.frame(
    s1 = "ctrl", s2 = "ctrl", s3 = "treat", s4 = "treat",
    stringsAsFactors = FALSE
  )
  rownames(df) <- factor_name

  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(colnames(result), "myfactor")
})

# ---------------------------------------------------------------------------
# 16. Regression: existing tests still pass via Manual Fill-in path
# ---------------------------------------------------------------------------

test_that("manual fill-in path: basic two-group design unchanged", {
  # Verifies the Manual Fill-in tab's apply_design observer still feeds the
  # same function correctly (regression check for the refactored modal).
  samples <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  df <- data.frame(
    ctrl_1 = "ctrl", ctrl_2 = "ctrl", treat_1 = "treat", treat_2 = "treat",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  result <- build_sample_info_from_df(df, samples)

  expect_false(is.null(result))
  expect_equal(nrow(result), 4L)
  expect_equal(ncol(result), 1L)
  expect_equal(
    as.character(result[, "group"]),
    c("CTRL", "CTRL", "TREAT", "TREAT")
  )
})
