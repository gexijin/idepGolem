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
  groups <- suppressWarnings(detect_groups(sample_names, max_length = 30))

  # Check that long group name was truncated to 30 chars
  long_group <- groups[1]
  expect_equal(nchar(long_group), 30)
  expect_equal(long_group, "VeryLongGroupNameThatExceeds30")

  # Check that short group name was not affected
  control_group <- unique(groups[3:4])
  expect_equal(control_group, "Control")
  expect_equal(nchar(control_group), 7)
})

test_that("detect_groups respects custom max_length parameter", {
  sample_names <- c(
    "VeryLongGroupName_1",
    "VeryLongGroupName_2",
    "Short_1",
    "Short_2"
  )

  # Test with max_length = 10
  groups <- suppressWarnings(detect_groups(sample_names, max_length = 10))

  long_group <- unique(groups[1:2])
  expect_equal(nchar(long_group), 10)
  expect_equal(long_group, "VeryLongGr")
})

test_that("detect_groups limits to 12 groups plus Other", {
  # Create 15 groups with varying frequencies
  sample_names <- c(
    paste0("GroupA_", 1:10),  # Most frequent
    paste0("GroupB_", 1:8),
    paste0("GroupC_", 1:6),
    paste0("GroupD_", 1:5),
    paste0("GroupE_", 1:4),
    paste0("GroupF_", 1:4),
    paste0("GroupG_", 1:3),
    paste0("GroupH_", 1:3),
    paste0("GroupI_", 1:3),
    paste0("GroupJ_", 1:2),
    paste0("GroupK_", 1:2),
    paste0("GroupL_", 1:2),
    paste0("GroupM_", 1:2),  # 12th most frequent
    paste0("GroupN_", 1:2),  # Should become Other
    paste0("GroupO_", 1:1)   # Should become Other
  )

  groups <- suppressWarnings(detect_groups(sample_names))

  # Should have at most 13 unique groups (12 + Other)
  expect_lte(length(unique(groups)), 13)

  # Should have "Other" group
  expect_true("Other" %in% groups)

  # Less frequent groups should be recoded as Other
  # GroupN and GroupO should not appear as separate groups
  expect_false("GroupN" %in% unique(groups))
  expect_false("GroupO" %in% unique(groups))
})

test_that("detect_groups preserves all groups when 12 or fewer", {
  # Create exactly 12 groups
  sample_names <- c(
    paste0("A_", 1:2), paste0("B_", 1:2), paste0("C_", 1:2),
    paste0("D_", 1:2), paste0("E_", 1:2), paste0("F_", 1:2),
    paste0("G_", 1:2), paste0("H_", 1:2), paste0("I_", 1:2),
    paste0("J_", 1:2), paste0("K_", 1:2), paste0("L_", 1:2)
  )

  groups <- suppressWarnings(detect_groups(sample_names))

  # Should have exactly 12 groups
  expect_equal(length(unique(groups)), 12)

  # Should NOT have "Other" group (all groups preserved)
  expect_false("Other" %in% groups)
})

test_that("detect_groups respects custom max_groups parameter", {
  # Create 8 groups, test with max_groups = 5
  sample_names <- c(
    paste0("A_", 1:3), paste0("B_", 1:3), paste0("C_", 1:3),
    paste0("D_", 1:3), paste0("E_", 1:3), paste0("F_", 1:3),
    paste0("G_", 1:3), paste0("H_", 1:3)
  )

  groups <- suppressWarnings(detect_groups(sample_names, max_groups = 5))

  # Should have at most 6 unique groups (5 + Other)
  expect_lte(length(unique(groups)), 6)

  # Should have "Other" group
  expect_true("Other" %in% groups)
})
