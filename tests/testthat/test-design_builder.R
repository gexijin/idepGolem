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

# --- near_duplicate_pairs() ---------------------------------------------------

test_that("near_duplicate_pairs flags distance-1 typos", {
  result <- near_duplicate_pairs(c("CTRL", "CTR", "TREATED"))
  expect_length(result, 1L)
  expect_true(grepl("CTRL", result) && grepl("CTR", result))
})

test_that("near_duplicate_pairs ignores exact duplicates and distance >= 2", {
  expect_length(near_duplicate_pairs(c("CTRL", "CTRL", "TREATED")), 0L)
  expect_length(near_duplicate_pairs(c("CTRL", "DRUG", "TREATED")), 0L)
})

test_that("near_duplicate_pairs is case-insensitive", {
  expect_length(near_duplicate_pairs(c("Ctrl", "ctr")), 1L)
})

test_that("near_duplicate_pairs handles edge inputs", {
  expect_length(near_duplicate_pairs(character(0)), 0L)
  expect_length(near_duplicate_pairs("ALONE"), 0L)
  expect_length(near_duplicate_pairs(c(NA_character_, "")), 0L)
  expect_length(near_duplicate_pairs(c("CTRL", NA_character_, "")), 0L)
})

test_that("near_duplicate_pairs honors max_dist argument", {
  # "kitten" -> "kitin" has distance 2
  expect_length(near_duplicate_pairs(c("kitten", "kitin"), max_dist = 1L), 0L)
  expect_length(near_duplicate_pairs(c("kitten", "kitin"), max_dist = 2L), 1L)
})

test_that("near_duplicate_pairs reports multiple pair hits", {
  # Distance-1 pairs: CTRL↔CTR, CTR↔TR, TRT↔TR
  # CTR↔TR and TRT↔TR are suppressed because max(nchar) <= 3 is OK but
  # actually max here is 3 which is > 2, so they're flagged. CTRL↔CTR has
  # max length 4, no shared digit-suffix root → also flagged. All 3 kept.
  result <- near_duplicate_pairs(c("CTRL", "CTR", "TRT", "TR"))
  expect_length(result, 3L)
})

test_that("near_duplicate_pairs skips numbered-series patterns", {
  # Common lab labeling — these are intentional, not typos.
  expect_length(near_duplicate_pairs(c("S1", "S2", "S3")), 0L)
  expect_length(near_duplicate_pairs(c("Rep1", "Rep2", "Rep3")), 0L)
  expect_length(near_duplicate_pairs(c("Sample1", "Sample10")), 0L)
  expect_length(near_duplicate_pairs(c("D1", "D2")), 0L)
  # Pure pair indices.
  expect_length(near_duplicate_pairs(c("1", "2", "3")), 0L)
})

test_that("near_duplicate_pairs skips very short labels (max length <= 2)", {
  # Length-2 strings: distance 1 = 50% mismatch, far past typo signal.
  expect_length(near_duplicate_pairs(c("WT", "KT")), 0L)
  expect_length(near_duplicate_pairs(c("A", "B")), 0L)
})

test_that("near_duplicate_pairs still flags genuine typos in longer labels", {
  expect_length(near_duplicate_pairs(c("CTRL", "CTR")), 1L)
  expect_length(near_duplicate_pairs(c("treated", "treatd")), 1L)
  # Different root, both length >= 3 → real typo signal.
  expect_length(near_duplicate_pairs(c("Day1", "Bay1")), 1L)
})

# --- has_non_ascii() ----------------------------------------------------------

test_that("has_non_ascii returns FALSE for pure-ASCII input", {
  expect_false(has_non_ascii(c("CTRL", "TREATED", "WT_REP1")))
  expect_false(has_non_ascii("hello123"))
})

test_that("has_non_ascii detects accented Latin, CJK, and emoji", {
  expect_true(has_non_ascii("café"))
  expect_true(has_non_ascii("北京"))
  expect_true(has_non_ascii("a🙂"))
  expect_true(has_non_ascii(c("CTRL", "café")))
})

test_that("has_non_ascii handles edge inputs", {
  expect_false(has_non_ascii(NULL))
  expect_false(has_non_ascii(character(0)))
  expect_false(has_non_ascii(NA_character_))
  expect_false(has_non_ascii(""))
  expect_false(has_non_ascii(c(NA_character_, "")))
})

# --- sanitize_created_names() (strict wizard sanitizer) ----------------------

test_that("sanitize_created_names keeps ASCII alphanumerics and underscores", {
  expect_equal(sanitize_created_names("WT_Rep1"), "WT_Rep1")
  expect_equal(sanitize_created_names("abc123"), "abc123")
  expect_equal(sanitize_created_names("foo_BAR_42"), "foo_BAR_42")
})

test_that("sanitize_created_names strips punctuation that sanitize_names preserves", {
  expect_equal(sanitize_created_names('what?'), "what")
  expect_equal(sanitize_created_names('say "hi"'), "sayhi")
  expect_equal(sanitize_created_names("ctrl(rep1)"), "ctrlrep1")
  expect_equal(sanitize_created_names("a!b@c#d$"), "abcd")
  expect_equal(sanitize_created_names("foo:bar"), "foobar")
})

test_that("sanitize_created_names strips spaces, hyphens, dots", {
  expect_equal(sanitize_created_names("Sample 1"), "Sample1")
  expect_equal(sanitize_created_names("Sample-1"), "Sample1")
  expect_equal(sanitize_created_names("Sample.1"), "Sample1")
})

test_that("sanitize_created_names strips non-ASCII characters", {
  expect_equal(sanitize_created_names("café"), "caf")
  expect_equal(sanitize_created_names("北京"), "")
  expect_equal(sanitize_created_names("a🙂b"), "ab")
})

test_that("sanitize_created_names preserves NA and handles edge inputs", {
  expect_equal(sanitize_created_names(NA_character_), NA_character_)
  expect_equal(
    sanitize_created_names(c("ok", NA_character_, "bad?")),
    c("ok", NA_character_, "bad")
  )
  expect_null(sanitize_created_names(NULL))
  expect_equal(sanitize_created_names(character(0)), character(0))
})

# --- sanitize_names() (original gentle behavior — unchanged contract) --------

test_that("sanitize_names preserves punctuation that the strict variant strips", {
  # The original sanitize_names only strips control chars, zero-width, NBSP,
  # hyphens, dots, and whitespace. Everything else passes through.
  expect_equal(sanitize_names('what?'), "what?")
  expect_equal(sanitize_names("ctrl(rep1)"), "ctrl(rep1)")
  expect_equal(sanitize_names("foo:bar"), "foo:bar")
})

test_that("sanitize_names still strips spaces, hyphens, dots (legacy behavior)", {
  expect_equal(sanitize_names("Sample 1"), "Sample1")
  expect_equal(sanitize_names("Sample-1"), "Sample1")
  expect_equal(sanitize_names("Sample.1"), "Sample1")
})

test_that("sanitize_names preserves non-ASCII characters", {
  # Uploaded-file sanitizer is intentionally permissive on Unicode.
  expect_equal(sanitize_names("café"), "café")
  expect_equal(sanitize_names("北京"), "北京")
})
