#!/usr/bin/env Rscript
# Measure production startup time for idepGolem

cat("\n")
cat("========================================\n")
cat("PRODUCTION STARTUP TIMING TEST\n")
cat("========================================\n")

.start_time <- Sys.time()
cat(sprintf("[%s] Loading idepGolem package...\n", format(.start_time, "%H:%M:%OS3")))

library(idepGolem)

.loaded_time <- Sys.time()
cat(sprintf("[%s] Package loaded (%.3fs)\n",
            format(.loaded_time, "%H:%M:%OS3"),
            as.numeric(difftime(.loaded_time, .start_time, units = "secs"))))

cat(sprintf("[%s] Starting run_app()...\n", format(Sys.time(), "%H:%M:%OS3")))

# Run the app
run_app()
