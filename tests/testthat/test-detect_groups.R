test_that("detect_groups identifies good grouping patterns", {
  sample_names <- c(
    "Control_1", "Control_2", "Control_3",
    "Treatment_1", "Treatment_2", "Treatment_3"
  )
  groups <- detect_groups(sample_names)

  expect_equal(length(unique(groups)), 2)
  expect_equal(unique(groups), c("Control", "Treatment"))
})

test_that("detect_groups handles no grouping pattern (each sample unique)", {
  sample_names <- c("SampleA", "SampleB", "SampleC", "SampleD", "SampleE")
  groups <- detect_groups(sample_names)

  # When 100% of samples are unique (>= 80% threshold), should return "Samples"
  expect_equal(length(unique(groups)), 1)
  expect_equal(unique(groups), "Samples")
})

test_that("detect_groups handles edge case at 80% threshold", {
  # 4 out of 5 groups are unique (80%)
  sample_names <- c("A1", "A2", "B", "C", "D")
  groups <- detect_groups(sample_names)

  # Should trigger threshold and return "Samples"
  expect_equal(length(unique(groups)), 1)
  expect_equal(unique(groups), "Samples")
})

test_that("detect_groups preserves groups below 80% threshold", {
  # 3 out of 6 groups are unique (50%)
  sample_names <- c("A1", "A2", "B1", "B2", "C1", "C2")
  groups <- detect_groups(sample_names)

  # Should keep the groups
  expect_equal(length(unique(groups)), 3)
  expect_true(all(unique(groups) %in% c("A", "B", "C")))
})

test_that("detect_groups works with sample_info matrix", {
  sample_names <- c("S1", "S2", "S3", "S4")
  sample_info <- data.frame(
    Group = c("Control", "Control", "Treatment", "Treatment"),
    row.names = sample_names
  )

  groups <- detect_groups(sample_names, sample_info)

  expect_equal(length(unique(groups)), 2)
  expect_true(all(unique(groups) %in% c("Control", "Treatment")))
})

test_that("detect_groups handles sample_info with unique combinations", {
  sample_names <- c("S1", "S2", "S3", "S4", "S5")
  # Create sample_info where each sample has unique factor combination
  sample_info <- data.frame(
    Factor1 = c("A", "B", "C", "D", "E"),
    Factor2 = c("X", "Y", "Z", "W", "V"),
    row.names = sample_names
  )

  groups <- detect_groups(sample_names, sample_info)

  # Should trigger threshold and return "Samples"
  expect_equal(length(unique(groups)), 1)
  expect_equal(unique(groups), "Samples")
})
