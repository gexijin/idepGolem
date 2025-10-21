#!/usr/bin/env Rscript
# Benchmark different CSV reading methods for demo_file_info.csv

cat("\n")
cat("========================================\n")
cat("CSV READING PERFORMANCE BENCHMARK\n")
cat("========================================\n")

# Set up the file path
db_ver <- "data113"
DATAPATH <- Sys.getenv("IDEP_DATABASE")[1]
if (nchar(DATAPATH) == 0) {
  DATAPATH <- paste0("../../data/")
}
DATAPATH <- paste0(DATAPATH, "/", db_ver, "/")

# If that doesn't exist, try local
if (!dir.exists(DATAPATH)) {
  DATAPATH <- paste0("./", db_ver, "/")
}

csv_file <- paste0(DATAPATH, "demo/demo_data_info.csv")

cat(sprintf("Testing file: %s\n", csv_file))

if (!file.exists(csv_file)) {
  cat("ERROR: File does not exist!\n")
  cat(sprintf("Checked path: %s\n", csv_file))
  quit(status = 1)
}

# Check file size
file_size <- file.info(csv_file)$size
cat(sprintf("File size: %d bytes (%.2f KB)\n\n", file_size, file_size/1024))

# Number of iterations for benchmarking
n_iterations <- 100

cat(sprintf("Running %d iterations of each method...\n\n", n_iterations))

# ============================================================================
# Method 1: base R read.csv() - Current method
# ============================================================================
cat("Method 1: read.csv() [CURRENT METHOD]\n")
cat("--------------------------------------\n")

times_read_csv <- numeric(n_iterations)
for (i in 1:n_iterations) {
  start <- Sys.time()
  data1 <- read.csv(csv_file)
  times_read_csv[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
}

cat(sprintf("  Mean:   %.6f seconds\n", mean(times_read_csv)))
cat(sprintf("  Median: %.6f seconds\n", median(times_read_csv)))
cat(sprintf("  Min:    %.6f seconds\n", min(times_read_csv)))
cat(sprintf("  Max:    %.6f seconds\n", max(times_read_csv)))
cat(sprintf("  Rows:   %d\n", nrow(data1)))
cat(sprintf("  Cols:   %d\n\n", ncol(data1)))

# ============================================================================
# Method 2: data.table::fread()
# ============================================================================
if (requireNamespace("data.table", quietly = TRUE)) {
  cat("Method 2: data.table::fread()\n")
  cat("--------------------------------------\n")

  times_fread <- numeric(n_iterations)
  for (i in 1:n_iterations) {
    start <- Sys.time()
    data2 <- data.table::fread(csv_file, data.table = FALSE)
    times_fread[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  cat(sprintf("  Mean:   %.6f seconds\n", mean(times_fread)))
  cat(sprintf("  Median: %.6f seconds\n", median(times_fread)))
  cat(sprintf("  Min:    %.6f seconds\n", min(times_fread)))
  cat(sprintf("  Max:    %.6f seconds\n", max(times_fread)))
  cat(sprintf("  Speedup: %.2fx faster than read.csv()\n\n", mean(times_read_csv) / mean(times_fread)))
} else {
  cat("Method 2: data.table::fread() - SKIPPED (package not installed)\n\n")
}

# ============================================================================
# Method 3: readr::read_csv()
# ============================================================================
if (requireNamespace("readr", quietly = TRUE)) {
  cat("Method 3: readr::read_csv()\n")
  cat("--------------------------------------\n")

  times_read_csv_readr <- numeric(n_iterations)
  for (i in 1:n_iterations) {
    start <- Sys.time()
    data3 <- readr::read_csv(csv_file, show_col_types = FALSE)
    times_read_csv_readr[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  cat(sprintf("  Mean:   %.6f seconds\n", mean(times_read_csv_readr)))
  cat(sprintf("  Median: %.6f seconds\n", median(times_read_csv_readr)))
  cat(sprintf("  Min:    %.6f seconds\n", min(times_read_csv_readr)))
  cat(sprintf("  Max:    %.6f seconds\n", max(times_read_csv_readr)))
  cat(sprintf("  Speedup: %.2fx faster than read.csv()\n\n", mean(times_read_csv) / mean(times_read_csv_readr)))
} else {
  cat("Method 3: readr::read_csv() - SKIPPED (package not installed)\n\n")
}

# ============================================================================
# Method 4: vroom::vroom()
# ============================================================================
if (requireNamespace("vroom", quietly = TRUE)) {
  cat("Method 4: vroom::vroom()\n")
  cat("--------------------------------------\n")

  times_vroom <- numeric(n_iterations)
  for (i in 1:n_iterations) {
    start <- Sys.time()
    data4 <- vroom::vroom(csv_file, show_col_types = FALSE)
    times_vroom[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  cat(sprintf("  Mean:   %.6f seconds\n", mean(times_vroom)))
  cat(sprintf("  Median: %.6f seconds\n", median(times_vroom)))
  cat(sprintf("  Min:    %.6f seconds\n", min(times_vroom)))
  cat(sprintf("  Max:    %.6f seconds\n", max(times_vroom)))
  cat(sprintf("  Speedup: %.2fx faster than read.csv()\n\n", mean(times_read_csv) / mean(times_vroom)))
} else {
  cat("Method 4: vroom::vroom() - SKIPPED (package not installed)\n\n")
}

# ============================================================================
# Summary Table
# ============================================================================
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

results <- data.frame(
  Method = "read.csv()",
  Mean_Time_s = mean(times_read_csv),
  Speedup = 1.0,
  stringsAsFactors = FALSE
)

if (exists("times_fread")) {
  results <- rbind(results, data.frame(
    Method = "fread()",
    Mean_Time_s = mean(times_fread),
    Speedup = mean(times_read_csv) / mean(times_fread),
    stringsAsFactors = FALSE
  ))
}

if (exists("times_read_csv_readr")) {
  results <- rbind(results, data.frame(
    Method = "read_csv()",
    Mean_Time_s = mean(times_read_csv_readr),
    Speedup = mean(times_read_csv) / mean(times_read_csv_readr),
    stringsAsFactors = FALSE
  ))
}

if (exists("times_vroom")) {
  results <- rbind(results, data.frame(
    Method = "vroom()",
    Mean_Time_s = mean(times_vroom),
    Speedup = mean(times_read_csv) / mean(times_vroom),
    stringsAsFactors = FALSE
  ))
}

# Sort by speed
results <- results[order(results$Mean_Time_s), ]

print(results, row.names = FALSE)

cat("\n")
cat("========================================\n")
cat("RECOMMENDATION\n")
cat("========================================\n")

fastest <- results[1, ]
cat(sprintf("Fastest method: %s\n", fastest$Method))
cat(sprintf("Time: %.6f seconds\n", fastest$Mean_Time_s))
cat(sprintf("Speedup: %.2fx faster than read.csv()\n", fastest$Speedup))

if (fastest$Method != "read.csv()") {
  cat(sprintf("\nRecommendation: Replace read.csv() with %s\n", fastest$Method))
  cat("This is a low-risk change that will improve startup performance.\n")
}

cat("\n")
