test_that("remove_ensembl_version handles single ID", {
  result <- remove_ensembl_version("ENSG00000211459.2")
  expect_equal(result, "ENSG00000211459")
})

test_that("remove_ensembl_version handles multiple IDs", {
  input <- c("ENSG00000211459.2", "ENSG00000123456.10", "ENSG00098765432.1")
  expected <- c("ENSG00000211459", "ENSG00000123456", "ENSG00098765432")
  result <- remove_ensembl_version(input)
  expect_equal(result, expected)
})

test_that("remove_ensembl_version handles IDs without versions", {
  input <- "ENSG00000211459"
  result <- remove_ensembl_version(input)
  expect_equal(result, input)
})

test_that("remove_ensembl_version handles mixed IDs", {
  input <- c("ENSG00000211459.2", "ENSG00000123456", "ENSG00000987654.15")
  expected <- c("ENSG00000211459", "ENSG00000123456", "ENSG00000987654")
  result <- remove_ensembl_version(input)
  expect_equal(result, expected)
})

test_that("remove_ensembl_version handles versions with more than 2 digits", {
  # Version with 3+ digits should NOT be removed
  result <- remove_ensembl_version("ENSG00000211459.100")
  expect_equal(result, "ENSG00000211459.100")
})

test_that("remove_ensembl_version handles NA values", {
  input <- c("ENSG00000211459.2", NA_character_, "ENSG00000123456.1")
  result <- remove_ensembl_version(input)
  expect_equal(result[1], "ENSG00000211459")
  expect_true(is.na(result[2]))
  expect_equal(result[3], "ENSG00000123456")
})

test_that("remove_ensembl_version handles empty strings", {
  result <- remove_ensembl_version("")
  expect_equal(result, "")
})

test_that("remove_ensembl_version handles different digit lengths", {
  input <- c(
    "ENSG00000001.5",      # 8 digits
    "ENSG0000000123.7",    # 10 digits
    "ENSG000000012345.2"   # 12 digits
  )
  expected <- c(
    "ENSG00000001",
    "ENSG0000000123",
    "ENSG000000012345"
  )
  result <- remove_ensembl_version(input)
  expect_equal(result, expected)
})

test_that("remove_ensembl_version is vectorized and efficient", {
  # Test with 10000 IDs to verify vectorization
  n <- 10000
  input <- paste0("ENSG", sprintf("%011d", 1:n), ".", sample(1:99, n, replace = TRUE))

  # Should complete quickly (well under 1 second for vectorized operation)
  start_time <- Sys.time()
  result <- remove_ensembl_version(input)
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  # Verify all versions were removed
  expect_true(all(!grepl("\\.", result)))

  # Verify length is preserved
  expect_equal(length(result), n)

  # Performance check: should process >10000 IDs per second (very conservative)
  expect_lt(elapsed, 1.0)
})

test_that("remove_ensembl_version leaves non-Ensembl strings unchanged", {
  # These strings should NOT be modified
  input <- c("ABC.cded", "GENE123.4", "XYZ.123", "random.text")
  result <- remove_ensembl_version(input)
  expect_equal(result, input)
})

test_that("remove_ensembl_version handles invalid Ensembl-like IDs", {
  # Too few digits (less than 8)
  expect_equal(remove_ensembl_version("ENSG1234567.2"), "ENSG1234567.2")
  expect_equal(remove_ensembl_version("ENSG123.4"), "ENSG123.4")

  # Too many digits (more than 12)
  expect_equal(remove_ensembl_version("ENSG1234567890123.5"), "ENSG1234567890123.5")

  # Not starting with ENSG
  expect_equal(remove_ensembl_version("ENST00000211459.2"), "ENST00000211459.2")
  expect_equal(remove_ensembl_version("ENS00000211459.2"), "ENS00000211459.2")
})

test_that("remove_ensembl_version handles mixed valid and invalid IDs", {
  input <- c(
    "ENSG00000211459.2",    # Valid Ensembl with version
    "ABC.cded",              # Non-Ensembl string
    "ENSG00000123456",       # Valid Ensembl without version
    "GENE123.4",             # Invalid (not ENSG)
    "ENSG123.5",             # Invalid (too few digits)
    "ENSG00000987654.15"     # Valid Ensembl with version
  )
  expected <- c(
    "ENSG00000211459",       # Version removed
    "ABC.cded",              # Unchanged
    "ENSG00000123456",       # Unchanged
    "GENE123.4",             # Unchanged
    "ENSG123.5",             # Unchanged
    "ENSG00000987654"        # Version removed
  )
  result <- remove_ensembl_version(input)
  expect_equal(result, expected)
})

test_that("remove_ensembl_version validates digit length boundaries", {
  # Exactly 8 digits (minimum) - should work
  expect_equal(remove_ensembl_version("ENSG12345678.2"), "ENSG12345678")

  # Exactly 12 digits (maximum) - should work
  expect_equal(remove_ensembl_version("ENSG123456789012.5"), "ENSG123456789012")

  # 7 digits - should NOT work
  expect_equal(remove_ensembl_version("ENSG1234567.2"), "ENSG1234567.2")

  # 13 digits - should NOT work
  expect_equal(remove_ensembl_version("ENSG1234567890123.2"), "ENSG1234567890123.2")
})

test_that("remove_ensembl_version handles dots without versions", {
  # String with dot but no digits after - dot should be removed
  input <- "ENSG00000211459."
  result <- remove_ensembl_version(input)
  expect_equal(result, "ENSG00000211459")  # Dot should be removed

  # Multiple IDs with trailing dots
  input2 <- c("ENSG00000111111.", "ENSG00000222222.", "ENSG00000333333")
  expected2 <- c("ENSG00000111111", "ENSG00000222222", "ENSG00000333333")
  result2 <- remove_ensembl_version(input2)
  expect_equal(result2, expected2)
})

test_that("remove_ensembl_version handles version numbers of varying length", {
  input <- c(
    "ENSG00000211459.1",     # 1 digit version - should be removed
    "ENSG00000211459.12",    # 2 digit version - should be removed
    "ENSG00000211459.123",   # 3 digit version - should NOT be removed
    "ENSG00000211459.1234"   # 4 digit version - should NOT be removed
  )
  expected <- c(
    "ENSG00000211459",       # 1-digit version removed
    "ENSG00000211459",       # 2-digit version removed
    "ENSG00000211459.123",   # 3-digit version kept
    "ENSG00000211459.1234"   # 4-digit version kept
  )
  result <- remove_ensembl_version(input)
  expect_equal(result, expected)
})

test_that("remove_ensembl_version handles whitespace correctly", {
  # IDs with whitespace should not match
  expect_equal(remove_ensembl_version(" ENSG00000211459.2"), " ENSG00000211459.2")
  expect_equal(remove_ensembl_version("ENSG00000211459.2 "), "ENSG00000211459.2 ")
  expect_equal(remove_ensembl_version("ENSG 00000211459.2"), "ENSG 00000211459.2")
})

test_that("remove_ensembl_version keeps 3-digit versions like .333", {
  # Specific test for user's example
  expect_equal(remove_ensembl_version("ENSG00000222222.333"), "ENSG00000222222.333")

  # More 3-digit examples
  expect_equal(remove_ensembl_version("ENSG00000111111.999"), "ENSG00000111111.999")
  expect_equal(remove_ensembl_version("ENSG00000555555.100"), "ENSG00000555555.100")
})

test_that("remove_ensembl_version performance with large mixed dataset", {
  # Test with 100000 mixed IDs (valid and invalid)
  n <- 100000
  valid_ids <- paste0("ENSG", sprintf("%011d", 1:(n/2)), ".", sample(1:99, n/2, replace = TRUE))
  invalid_ids <- c(
    rep("ABC.cded", n/4),
    rep("GENE123.4", n/4)
  )
  input <- c(valid_ids, invalid_ids)

  # Should still be very fast
  start_time <- Sys.time()
  result <- remove_ensembl_version(input)
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  # Verify length is preserved
  expect_equal(length(result), n)

  # Verify valid IDs had versions removed
  expect_true(all(!grepl("\\.", result[1:(n/2)])))

  # Verify invalid IDs were unchanged
  expect_equal(result[(n/2 + 1):n], invalid_ids)

  # Performance check: should process >50000 IDs per second (conservative)
  expect_lt(elapsed, 2.0)
})
