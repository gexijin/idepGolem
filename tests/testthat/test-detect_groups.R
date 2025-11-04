test_that("detect_groups identifies good grouping patterns", {
  sample_names <- c(
    "Control_1", "Control_2", "Control_3",
    "Treatment_1", "Treatment_2", "Treatment_3"
  )
  groups <- detect_groups(sample_names)

  expect_equal(length(unique(groups)), 2)
  expect_equal(sort(unique(groups)), c("Control", "Treatment"))
})

test_that("detect_groups handles no grouping pattern (all samples unique)", {
  sample_names <- c("SampleA", "SampleB", "SampleC", "SampleD", "SampleE")
  groups <- detect_groups(sample_names)

  # When all samples are unique, all become "Other", then all become "Samples"
  expect_equal(length(unique(groups)), 1)
  expect_equal(unique(groups), "Samples")
})

test_that("detect_groups labels non-replicated groups as Other", {
  # Mix of replicated and non-replicated groups
  sample_names <- c("A1", "A2", "B1", "B2", "C", "D", "E")
  groups <- detect_groups(sample_names)

  # A and B have replicates, C/D/E are unique -> become "Other"
  expect_equal(length(unique(groups)), 3)
  expect_true(all(unique(groups) %in% c("A", "B", "Other")))
  expect_equal(sum(groups == "Other"), 3) # C, D, E
})

test_that("detect_groups preserves all groups when all have replicates", {
  sample_names <- c("A1", "A2", "B1", "B2", "C1", "C2")
  groups <- detect_groups(sample_names)

  # All groups have replicates, should keep all
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

test_that("detect_groups truncates long group names", {
  # Create samples with group names longer than 30 characters
  sample_names <- c(
    "VeryLongGroupNameThatExceeds30Characters_1",
    "VeryLongGroupNameThatExceeds30Characters_2",
    "Control_1",
    "Control_2"
  )

  # Suppress warnings for testing
  groups <- suppressWarnings(detect_groups(sample_names))

  # Check that long group name was truncated to 30 chars
  long_group <- groups[1]
  expect_equal(nchar(long_group), 30)
  expect_equal(long_group, "VeryLongGroupNameThatExceeds30")

  # Check that short group name was not affected
  control_group <- unique(groups[3:4])
  expect_equal(control_group, "Control")
  expect_equal(nchar(control_group), 7)
})
