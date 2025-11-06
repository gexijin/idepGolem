test_that("remove_gene_version handles single ID", {
  result <- remove_gene_version("ENSG00000211459.2")
  expect_equal(result, "ENSG00000211459")
})

test_that("remove_gene_version handles multiple IDs", {
  input <- c("ENSG00000211459.2", "ENSG00000123456.10", "ENSG00098765432.1")
  expected <- c("ENSG00000211459", "ENSG00000123456", "ENSG00098765432")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles IDs without versions", {
  input <- "ENSG00000211459"
  result <- remove_gene_version(input)
  expect_equal(result, input)
})

test_that("remove_gene_version handles mixed IDs", {
  input <- c("ENSG00000211459.2", "ENSG00000123456", "ENSG00000987654.15")
  expected <- c("ENSG00000211459", "ENSG00000123456", "ENSG00000987654")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles versions with more than 2 digits", {
  # Version with 3+ digits should NOT be removed
  result <- remove_gene_version("ENSG00000211459.100")
  expect_equal(result, "ENSG00000211459.100")
})

test_that("remove_gene_version handles NA values", {
  input <- c("ENSG00000211459.2", NA_character_, "ENSG00000123456.1")
  result <- remove_gene_version(input)
  expect_equal(result[1], "ENSG00000211459")
  expect_true(is.na(result[2]))
  expect_equal(result[3], "ENSG00000123456")
})

test_that("remove_gene_version handles empty strings", {
  result <- remove_gene_version("")
  expect_equal(result, "")
})

test_that("remove_gene_version requires exactly 11 digits for Ensembl", {
  # Only IDs with exactly 11 digits should have versions removed
  input <- c(
    "ENSG00000000001.5",   # 11 digits - should work
    "ENSG0000000001.5",    # 10 digits - should NOT work
    "ENSG000000000001.2"   # 12 digits - should NOT work
  )
  expected <- c(
    "ENSG00000000001",           # Version removed (11 digits)
    "ENSG0000000001.5",          # Unchanged (10 digits)
    "ENSG000000000001.2"         # Unchanged (12 digits)
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version is vectorized and efficient", {
  # Test with 10000 IDs to verify vectorization
  n <- 10000
  input <- paste0("ENSG", sprintf("%011d", 1:n), ".", sample(1:99, n, replace = TRUE))

  # Should complete quickly (well under 1 second for vectorized operation)
  start_time <- Sys.time()
  result <- remove_gene_version(input)
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  # Verify all versions were removed
  expect_true(all(!grepl("\\.", result)))

  # Verify length is preserved
  expect_equal(length(result), n)

  # Performance check: should process >10000 IDs per second (very conservative)
  expect_lt(elapsed, 1.0)
})

test_that("remove_gene_version leaves non-Ensembl strings unchanged", {
  # These strings should NOT be modified
  input <- c("ABC.cded", "GENE123.4", "XYZ.123", "random.text")
  result <- remove_gene_version(input)
  expect_equal(result, input)
})

test_that("remove_gene_version handles invalid Ensembl-like IDs", {
  # Wrong digit count (not 11)
  expect_equal(remove_gene_version("ENSG1234567.2"), "ENSG1234567.2")
  expect_equal(remove_gene_version("ENSG123.4"), "ENSG123.4")
  expect_equal(remove_gene_version("ENSG1234567890123.5"), "ENSG1234567890123.5")

  # Missing letters after ENS (no species/type code)
  expect_equal(remove_gene_version("ENS00000211459.2"), "ENS00000211459.2")

  # Any letter combination is valid, just needs 11 digits
  # This should work since it has ENS + letter + 11 digits
  expect_equal(remove_gene_version("ENSX00000211459.2"), "ENSX00000211459")
})

test_that("remove_gene_version handles mixed valid and invalid IDs", {
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
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version validates exactly 11 digits for Ensembl", {
  # Exactly 11 digits - should work
  expect_equal(remove_gene_version("ENSG12345678901.2"), "ENSG12345678901")

  # 10 digits - should NOT work
  expect_equal(remove_gene_version("ENSG1234567890.2"), "ENSG1234567890.2")

  # 12 digits - should NOT work
  expect_equal(remove_gene_version("ENSG123456789012.2"), "ENSG123456789012.2")
})

test_that("remove_gene_version handles dots without versions", {
  # String with dot but no digits after - dot should be removed
  input <- "ENSG00000211459."
  result <- remove_gene_version(input)
  expect_equal(result, "ENSG00000211459")  # Dot should be removed

  # Multiple IDs with trailing dots
  input2 <- c("ENSG00000111111.", "ENSG00000222222.", "ENSG00000333333")
  expected2 <- c("ENSG00000111111", "ENSG00000222222", "ENSG00000333333")
  result2 <- remove_gene_version(input2)
  expect_equal(result2, expected2)
})

test_that("remove_gene_version handles version numbers of varying length", {
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
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles whitespace correctly", {
  # IDs with whitespace should not match
  expect_equal(remove_gene_version(" ENSG00000211459.2"), " ENSG00000211459.2")
  expect_equal(remove_gene_version("ENSG00000211459.2 "), "ENSG00000211459.2 ")
  expect_equal(remove_gene_version("ENSG 00000211459.2"), "ENSG 00000211459.2")
})

test_that("remove_gene_version keeps 3-digit versions like .333", {
  # Specific test for user's example
  expect_equal(remove_gene_version("ENSG00000222222.333"), "ENSG00000222222.333")

  # More 3-digit examples
  expect_equal(remove_gene_version("ENSG00000111111.999"), "ENSG00000111111.999")
  expect_equal(remove_gene_version("ENSG00000555555.100"), "ENSG00000555555.100")
})

test_that("remove_gene_version performance with large mixed dataset", {
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
  result <- remove_gene_version(input)
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

# Tests for RefSeq NM_ pattern support

test_that("remove_gene_version handles single RefSeq NM_ ID", {
  result <- remove_gene_version("NM_000546.5")
  expect_equal(result, "NM_000546")
})

test_that("remove_gene_version handles multiple NM_ IDs", {
  input <- c("NM_000546.5", "NM_001234.12", "NM_123456789.1")
  expected <- c("NM_000546", "NM_001234", "NM_123456789")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles mixed Ensembl and RefSeq IDs", {
  input <- c("ENSG00000211459.2", "NM_000546.5", "ENSG00000123456.10", "NM_001234.7")
  expected <- c("ENSG00000211459", "NM_000546", "ENSG00000123456", "NM_001234")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles NM_ IDs without versions", {
  input <- c("NM_000546", "NM_001234")
  result <- remove_gene_version(input)
  expect_equal(result, input)
})

test_that("remove_gene_version handles NM_ IDs with 3+ digit versions", {
  # Version with 3+ digits should NOT be removed
  input <- c("NM_000546.100", "NM_001234.333")
  result <- remove_gene_version(input)
  expect_equal(result, input)
})

test_that("remove_gene_version handles NM_ IDs with trailing dots", {
  input <- c("NM_000546.", "NM_001234.")
  expected <- c("NM_000546", "NM_001234")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version validates NM_ digit length boundaries", {
  # Exactly 6 digits (minimum) - should work
  expect_equal(remove_gene_version("NM_123456.2"), "NM_123456")

  # Exactly 9 digits (maximum) - should work
  expect_equal(remove_gene_version("NM_123456789.5"), "NM_123456789")

  # 5 digits - should NOT work
  expect_equal(remove_gene_version("NM_12345.2"), "NM_12345.2")

  # 10 digits - should NOT work
  expect_equal(remove_gene_version("NM_1234567890.2"), "NM_1234567890.2")
})

test_that("remove_gene_version handles case-insensitive NM_ matching", {
  input <- c("NM_000546.5", "nm_001234.7", "Nm_123456.2")
  expected <- c("NM_000546", "nm_001234", "Nm_123456")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version leaves invalid NM_-like IDs unchanged", {
  # Missing underscore
  expect_equal(remove_gene_version("NM000546.2"), "NM000546.2")

  # Unknown RefSeq prefix
  expect_equal(remove_gene_version("NT_000546.2"), "NT_000546.2")
  expect_equal(remove_gene_version("NG_000546.2"), "NG_000546.2")
})

test_that("remove_gene_version handles mixed valid Ensembl, RefSeq, and invalid IDs", {
  input <- c(
    "ENSG00000211459.2",    # Valid Ensembl
    "NM_000546.5",          # Valid RefSeq
    "ABC.cded",              # Invalid
    "NM_001234",             # Valid RefSeq without version
    "ENSG00000123456.10",    # Valid Ensembl
    "GENE123.4",             # Invalid
    "NM_12345.2",            # Invalid (too few digits)
    "ENSG123.5"              # Invalid (too few digits)
  )
  expected <- c(
    "ENSG00000211459",       # Version removed
    "NM_000546",             # Version removed
    "ABC.cded",              # Unchanged
    "NM_001234",             # Unchanged
    "ENSG00000123456",       # Version removed
    "GENE123.4",             # Unchanged
    "NM_12345.2",            # Unchanged
    "ENSG123.5"              # Unchanged
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles NM_ IDs with varying version lengths", {
  input <- c(
    "NM_000546.1",     # 1 digit version - should be removed
    "NM_000546.12",    # 2 digit version - should be removed
    "NM_000546.123",   # 3 digit version - should NOT be removed
    "NM_000546.1234"   # 4 digit version - should NOT be removed
  )
  expected <- c(
    "NM_000546",       # 1-digit version removed
    "NM_000546",       # 2-digit version removed
    "NM_000546.123",   # 3-digit version kept
    "NM_000546.1234"   # 4-digit version kept
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

# Tests for Ensembl transcript IDs (ENST)

test_that("remove_gene_version handles Ensembl transcript IDs", {
  input <- c("ENST00000456328.2", "ENST00000450305.5", "ENSMUST00000082423.1")
  expected <- c("ENST00000456328", "ENST00000450305", "ENSMUST00000082423")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles ENST without versions", {
  input <- "ENST00000456328"
  result <- remove_gene_version(input)
  expect_equal(result, input)
})

# Tests for Ensembl protein IDs (ENSP)

test_that("remove_gene_version handles Ensembl protein IDs", {
  input <- c("ENSP00000384458.1", "ENSP00000362111.3", "ENSMUSP00000082423.2")
  expected <- c("ENSP00000384458", "ENSP00000362111", "ENSMUSP00000082423")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

# Tests for species-specific Ensembl IDs

test_that("remove_gene_version handles species-specific Ensembl IDs", {
  input <- c(
    "ENSMUSG00000025902.5",   # Mouse gene
    "ENSMUST00000027073.9",   # Mouse transcript
    "ENSMUSP00000027073.4",   # Mouse protein
    "ENSDARG00000000001.10",  # Zebrafish gene
    "ENSRNOT00000000001.5"    # Rat transcript
  )
  expected <- c(
    "ENSMUSG00000025902",
    "ENSMUST00000027073",
    "ENSMUSP00000027073",
    "ENSDARG00000000001",
    "ENSRNOT00000000001"
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

# Tests for additional RefSeq ID types

test_that("remove_gene_version handles RefSeq protein IDs (NP_)", {
  input <- c("NP_000537.3", "NP_001234.5", "NP_123456789.12")
  expected <- c("NP_000537", "NP_001234", "NP_123456789")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles RefSeq ncRNA IDs (NR_)", {
  input <- c("NR_046018.2", "NR_001234.7", "NR_123456.1")
  expected <- c("NR_046018", "NR_001234", "NR_123456")
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

test_that("remove_gene_version handles predicted RefSeq IDs (XM_, XR_, XP_)", {
  input <- c(
    "XM_011545467.2",   # Predicted mRNA
    "XM_001234567.5",
    "XR_007058843.1",   # Predicted ncRNA
    "XR_001234567.3",
    "XP_011545467.1",   # Predicted protein
    "XP_001234567.3"
  )
  expected <- c(
    "XM_011545467",
    "XM_001234567",
    "XR_007058843",
    "XR_001234567",
    "XP_011545467",
    "XP_001234567"
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})

# Comprehensive mixed test

test_that("remove_gene_version handles all ID types mixed together", {
  input <- c(
    "ENSG00000211459.2",      # Human gene
    "ENSMUSG00000025902.5",   # Mouse gene
    "ENST00000456328.2",      # Transcript
    "ENSP00000384458.1",      # Protein
    "NM_000546.5",            # RefSeq mRNA
    "NP_000537.3",            # RefSeq protein
    "NR_046018.2",            # RefSeq ncRNA
    "XM_011545467.2",         # Predicted mRNA
    "XR_007058843.1",         # Predicted ncRNA
    "XP_011545467.1",         # Predicted protein
    "INVALID.123",            # Invalid ID
    "ENSG00000123456"         # No version
  )
  expected <- c(
    "ENSG00000211459",
    "ENSMUSG00000025902",
    "ENST00000456328",
    "ENSP00000384458",
    "NM_000546",
    "NP_000537",
    "NR_046018",
    "XM_011545467",
    "XR_007058843",
    "XP_011545467",
    "INVALID.123",
    "ENSG00000123456"
  )
  result <- remove_gene_version(input)
  expect_equal(result, expected)
})
