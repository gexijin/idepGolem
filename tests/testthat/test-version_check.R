test_that("version comparison logic works correctly", {
  # Mock function to test version comparison logic
  compare_versions <- function(current_ver, latest_tag) {
    current_parts <- as.numeric(strsplit(current_ver, "\\.")[[1]])
    latest_parts <- as.numeric(strsplit(latest_tag, "\\.")[[1]])

    if (length(current_parts) < 2 || length(latest_parts) < 2) {
      return(FALSE)
    }

    major_update <- latest_parts[1] > current_parts[1]
    minor_update <- latest_parts[1] == current_parts[1] && latest_parts[2] > current_parts[2]

    return(major_update || minor_update)
  }

  # Test major version update
  expect_true(compare_versions("2.20", "3.0"))
  expect_true(compare_versions("2.20.1", "3.0.0"))

  # Test minor version update
  expect_true(compare_versions("2.20", "2.21"))
  expect_true(compare_versions("2.20.3", "2.21.0"))

  # Test patch version update (should NOT trigger)
  expect_false(compare_versions("2.20.1", "2.20.2"))
  expect_false(compare_versions("2.20.3", "2.20.4"))

  # Test same version
  expect_false(compare_versions("2.20", "2.20"))
  expect_false(compare_versions("2.20.1", "2.20.1"))

  # Test older version
  expect_false(compare_versions("2.21", "2.20"))
  expect_false(compare_versions("3.0", "2.20"))
})
