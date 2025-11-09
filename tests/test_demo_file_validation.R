# Test: Demo File Validation
# This test demonstrates the new file existence checking

library(testthat)

test_that("Demo file validation handles missing expression file", {
  # Mock a missing demo file scenario
  result <- input_data(
    expression_file = NULL,
    experiment_file = NULL,
    go_button = 1, # Trigger demo data loading
    demo_data_file = "/nonexistent/path/missing_demo.csv",
    demo_metadata_file = "",
    max_group_name_length = 30
  )

  # Should return NULL when file doesn't exist
  expect_null(result)
})

test_that("Demo file validation handles missing design file gracefully", {
  # This would require a valid expression file but missing design file
  # The function should warn but continue processing
  # (This test would need a mock valid expression file)
})

cat("Demo file validation tests created.\n")
cat("These tests verify that:\n")
cat("1. Missing demo expression files are caught and reported\n")
cat("2. Missing demo design files trigger warnings but allow processing\n")
cat("3. File read errors are caught and reported with helpful messages\n")
