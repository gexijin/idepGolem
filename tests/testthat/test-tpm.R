test_that("compute_tpm produces column sums of 1e6", {
  counts <- matrix(
    c(10, 20,
      40, 80,
      0,  100),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))
  )
  lengths <- c(g1 = 1000, g2 = 2000, g3 = 500)

  result <- compute_tpm(counts, lengths)

  expect_equal(unname(colSums(result$tpm)), c(1e6, 1e6))
  expect_equal(result$n_dropped, 0L)
  expect_equal(rownames(result$tpm), c("g1", "g2", "g3"))
})

test_that("compute_tpm matches a hand-computed example", {
  # Two genes, equal counts of 100 in one sample.
  # g1 length 1 kb, g2 length 4 kb.
  # rpk = c(100, 25); scaling = 125; tpm = c(800000, 200000).
  counts <- matrix(
    c(100, 100),
    nrow = 2,
    dimnames = list(c("g1", "g2"), "s1")
  )
  lengths <- c(g1 = 1000, g2 = 4000)

  result <- compute_tpm(counts, lengths)

  expect_equal(unname(result$tpm[, 1]), c(800000, 200000))
})

test_that("compute_tpm drops genes with NA or zero length", {
  counts <- matrix(
    c(10, 20, 30, 40),
    nrow = 4, ncol = 1,
    dimnames = list(c("g1", "g2", "g3", "g4"), "s1")
  )
  lengths <- c(g1 = 1000, g2 = NA_real_, g3 = 0, g4 = 2000)

  result <- compute_tpm(counts, lengths)

  expect_equal(rownames(result$tpm), c("g1", "g4"))
  expect_equal(result$n_dropped, 2L)
  expect_equal(unname(colSums(result$tpm)), 1e6)
})

test_that("compute_tpm drops genes missing from lengths vector", {
  counts <- matrix(
    c(10, 20, 30),
    nrow = 3, ncol = 1,
    dimnames = list(c("g1", "g2", "g3"), "s1")
  )
  lengths <- c(g1 = 1000, g3 = 500)

  result <- compute_tpm(counts, lengths)

  expect_equal(sort(rownames(result$tpm)), c("g1", "g3"))
  expect_equal(result$n_dropped, 1L)
})

test_that("compute_tpm handles no overlap without crashing", {
  counts <- matrix(
    c(10, 20),
    nrow = 2, ncol = 1,
    dimnames = list(c("g1", "g2"), "s1")
  )
  lengths <- c(other_gene = 1000)

  result <- compute_tpm(counts, lengths)

  expect_equal(nrow(result$tpm), 0L)
  expect_equal(result$n_dropped, 2L)
})

test_that("compute_tpm handles empty counts matrix", {
  counts <- matrix(numeric(0), nrow = 0, ncol = 0)
  result <- compute_tpm(counts, c(g1 = 1000))

  expect_equal(nrow(result$tpm), 0L)
  expect_equal(result$n_dropped, 0L)
})
