test_that("input_data truncates long factor names when safe", {
  expr_file <- tempfile(fileext = ".csv")
  design_file <- tempfile(fileext = ".csv")

  expr_df <- data.frame(
    Gene = c("gene1", "gene2"),
    SampleA = c(10, 20),
    check.names = FALSE
  )

  design_df <- data.frame(
    SampleA = c("control", "control"),
    row.names = c(
      "ThisIsAVeryLongFactorNameThatNeedsTruncation",
      "ShortFactor"
    ),
    check.names = FALSE
  )

  write.csv(expr_df, expr_file, row.names = FALSE, quote = TRUE)
  write.csv(design_df, design_file, row.names = TRUE, quote = TRUE)

  result <- input_data(
    expression_file = list(datapath = expr_file),
    experiment_file = list(datapath = design_file),
    go_button = 0,
    demo_data_file = NULL,
    demo_metadata_file = NULL,
    max_group_name_length = 10
  )

  expect_false(is.null(result$sample_info))
  expect_equal(colnames(result$sample_info), c("ThisIsAVer", "ShortFactor"))
})

test_that("input_data keeps original factor names when truncation would duplicate", {
  expr_file <- tempfile(fileext = ".csv")
  design_file <- tempfile(fileext = ".csv")

  expr_df <- data.frame(
    Gene = c("gene1", "gene2"),
    SampleA = c(10, 20),
    check.names = FALSE
  )

  design_df <- data.frame(
    SampleA = c("control", "control"),
    row.names = c(
      "RepeatedFactorName_AAA",
      "RepeatedFactorName_BBB"
    ),
    check.names = FALSE
  )

  write.csv(expr_df, expr_file, row.names = FALSE, quote = TRUE)
  write.csv(design_df, design_file, row.names = TRUE, quote = TRUE)

  result <- input_data(
    expression_file = list(datapath = expr_file),
    experiment_file = list(datapath = design_file),
    go_button = 0,
    demo_data_file = NULL,
    demo_metadata_file = NULL,
    max_group_name_length = 10
  )

  expect_false(is.null(result$sample_info))
  expect_equal(
    colnames(result$sample_info),
    c("RepeatedFactorName_AAA", "RepeatedFactorName_BBB")
  )
})

test_that("input_data truncates factor levels when safe", {
  expr_file <- tempfile(fileext = ".csv")
  design_file <- tempfile(fileext = ".csv")

  expr_df <- data.frame(
    Gene = c("gene1", "gene2"),
    SampleA = c(10, 20),
    SampleB = c(15, 25),
    check.names = FALSE
  )

  design_df <- data.frame(
    SampleA = c("AlphaLevelNameEXTREMELYLONG", "BetaLevelNameEXTREMELYLONG"),
    SampleB = c("GammaLevelNameEXTREMELYLONG", "DeltaLevelNameEXTREMELYLONG"),
    row.names = c("FactorOne", "FactorTwo"),
    check.names = FALSE
  )

  write.csv(expr_df, expr_file, row.names = FALSE, quote = TRUE)
  write.csv(design_df, design_file, row.names = TRUE, quote = TRUE)

  result <- input_data(
    expression_file = list(datapath = expr_file),
    experiment_file = list(datapath = design_file),
    go_button = 0,
    demo_data_file = NULL,
    demo_metadata_file = NULL,
    max_group_name_length = 10
  )

  expect_false(is.null(result$sample_info))
  expect_equal(
    unname(result$sample_info[, "FactorOne"]),
    c("AlphaLevel", "GammaLevel")
  )
  expect_equal(
    unname(result$sample_info[, "FactorTwo"]),
    c("BetaLevelN", "DeltaLevel")
  )
})

test_that("input_data keeps factor levels when truncation would duplicate", {
  expr_file <- tempfile(fileext = ".csv")
  design_file <- tempfile(fileext = ".csv")

  expr_df <- data.frame(
    Gene = c("gene1", "gene2"),
    SampleA = c(10, 20),
    SampleB = c(15, 25),
    check.names = FALSE
  )

  design_df <- data.frame(
    SampleA = c("RepeatedLevelName_AAA", "UniqueLevel"),
    SampleB = c("RepeatedLevelName_BBB", "UniqueLevel"),
    row.names = c("FactorOne", "FactorTwo"),
    check.names = FALSE
  )

  write.csv(expr_df, expr_file, row.names = FALSE, quote = TRUE)
  write.csv(design_df, design_file, row.names = TRUE, quote = TRUE)

  result <- input_data(
    expression_file = list(datapath = expr_file),
    experiment_file = list(datapath = design_file),
    go_button = 0,
    demo_data_file = NULL,
    demo_metadata_file = NULL,
    max_group_name_length = 10
  )

  expect_false(is.null(result$sample_info))
  expect_equal(
    unname(result$sample_info[, "FactorOne"]),
    c("RepeatedLevelName_AAA", "RepeatedLevelName_BBB")
  )
})
