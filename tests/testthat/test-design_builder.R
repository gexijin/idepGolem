# Tests for create_design_template() and build_sample_info_from_df()
# These functions support the interactive design builder GUI (issue #363).

# --- create_design_template() -------------------------------------------------

test_that("create_design_template returns a 1-row data frame with correct structure", {
  sample_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
  result <- create_design_template(sample_names)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1L)
  expect_equal(colnames(result), sample_names)
  expect_equal(rownames(result), "group")
})

test_that("create_design_template pre-fills groups from descriptive sample names", {
  sample_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
  result <- create_design_template(sample_names)

  groups <- as.character(result[1, ])
  expect_equal(length(unique(groups)), 2L)
  expect_true(all(unique(groups) %in% c("Control", "Treatment")))
})

test_that("create_design_template handles generic sample names gracefully", {
  sample_names <- c("SRR001", "SRR002", "SRR003", "SRR004")
  result <- create_design_template(sample_names)

  expect_equal(nrow(result), 1L)
  expect_equal(colnames(result), sample_names)
  # All same group — detect_groups() falls back to "Samples"
  expect_equal(length(unique(as.character(result[1, ]))), 1L)
})

# --- build_sample_info_from_df() ----------------------------------------------

test_that("build_sample_info_from_df returns correct transposed matrix for single factor", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = "ctrl", s2 = "ctrl", s3 = "treat", s4 = "treat",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4L)   # rows = samples
  expect_equal(ncol(result), 1L)   # cols = factors
  expect_equal(rownames(result), sample_names)
  expect_equal(colnames(result), "group")
})

test_that("build_sample_info_from_df handles two-factor design", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = c("ctrl", "T0"), s2 = c("ctrl", "T1"),
    s3 = c("treat", "T0"), s4 = c("treat", "T1"),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c("group", "timepoint")

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  expect_equal(ncol(result), 2L)
  expect_equal(colnames(result), c("group", "timepoint"))
})

test_that("build_sample_info_from_df returns NULL when all samples share the same group", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = "ctrl", s2 = "ctrl", s3 = "ctrl", s4 = "ctrl",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  expect_null(build_sample_info_from_df(df, sample_names))
})

test_that("build_sample_info_from_df drops factors with no variation, keeps valid ones", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = c("ctrl", "same"), s2 = c("ctrl", "same"),
    s3 = c("treat", "same"), s4 = c("treat", "same"),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c("group", "boring")

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "group")
})

test_that("build_sample_info_from_df returns NULL when all factors have no variation", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = c("x", "y"), s2 = c("x", "y"), s3 = c("x", "y"), s4 = c("x", "y"),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c("f1", "f2")

  expect_null(build_sample_info_from_df(df, sample_names))
})

test_that("build_sample_info_from_df drops rows with empty factor names", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = c("ctrl", "A"), s2 = c("ctrl", "B"),
    s3 = c("treat", "A"), s4 = c("treat", "B"),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c("group", "")

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "group")
})

test_that("build_sample_info_from_df returns NULL for mismatched sample names", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    wrong1 = "ctrl", wrong2 = "ctrl", wrong3 = "treat", wrong4 = "treat",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  expect_null(build_sample_info_from_df(df, sample_names))
})

test_that("build_sample_info_from_df uppercases and strips spaces from cell values", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = "wild type", s2 = "wild type", s3 = "knock-out", s4 = "knock-out",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  cell_vals <- unique(as.character(result[, "group"]))
  expect_true(all(cell_vals %in% c("WILDTYPE", "KNOCKOUT")))
})

test_that("build_sample_info_from_df reorders columns to match sample_names order", {
  sample_names <- c("s1", "s2", "s3", "s4")
  # supply columns in reversed order
  df <- data.frame(
    s4 = "treat", s3 = "treat", s2 = "ctrl", s1 = "ctrl",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  result <- build_sample_info_from_df(df, sample_names)

  expect_false(is.null(result))
  expect_equal(rownames(result), sample_names)
  expect_equal(as.character(result[, "group"]), c("CTRL", "CTRL", "TREAT", "TREAT"))
})

test_that("build_sample_info_from_df truncates long labels when safe", {
  sample_names <- c("s1", "s2", "s3", "s4")
  df <- data.frame(
    s1 = "VeryLongGroupLabelHere", s2 = "VeryLongGroupLabelHere",
    s3 = "OtherLongGroupLabelHere", s4 = "OtherLongGroupLabelHere",
    stringsAsFactors = FALSE
  )
  rownames(df) <- "group"

  result <- build_sample_info_from_df(df, sample_names, max_group_name_length = 10)

  expect_false(is.null(result))
  expect_true(all(nchar(as.character(result[, "group"])) <= 10))
})

test_that("build_sample_info_from_df returns NULL for empty input", {
  expect_null(build_sample_info_from_df(NULL, c("s1", "s2")))
  expect_null(build_sample_info_from_df(
    data.frame(), c("s1", "s2")
  ))
})
