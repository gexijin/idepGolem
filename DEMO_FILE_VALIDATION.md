# Demo File Validation - Implementation Summary

## Problem
When a user selected a demo file that didn't exist on the server, the application would crash with an unhandled R error, providing a poor user experience.

## Solution
Added comprehensive file existence checking and error handling in `R/fct_01_load_data.R`.

## Changes Made

### 1. Demo Expression File Validation (Lines 150-167)
**Location:** `input_data()` function, after demo file path is set

**What it does:**
- Checks if the demo expression file exists before attempting to read it
- Shows a user-friendly error notification if file is missing
- Returns `NULL` gracefully instead of crashing
- Falls back to a warning message in non-Shiny contexts

**Code added:**
```r
# Validate demo file exists
if (!file.exists(in_file_data)) {
  if (requireNamespace("shiny", quietly = TRUE) && !is.null(shiny::getDefaultReactiveDomain())) {
    shiny::showNotification(
      ui = paste(
        "Error: Demo expression file not found:",
        basename(in_file_data),
        "Please contact the administrator or try a different demo dataset."
      ),
      id = "demo_file_not_found",
      duration = NULL,
      type = "error"
    )
  } else {
    warning("Demo expression file not found: ", in_file_data)
  }
  return(NULL)
}
```

### 2. Demo Design File Validation (Lines 322-347)
**Location:** `input_data()` function, when processing demo design files

**What it does:**
- Checks if the demo design file exists before attempting to read it
- Shows a warning (not error) if design file is missing
- **Allows processing to continue** without the design file (graceful degradation)
- Only processes file truncation logic if the design file loaded successfully

**Code added:**
```r
# Validate demo design file exists
if (!file.exists(demo_metadata_file)) {
  if (requireNamespace("shiny", quietly = TRUE) && !is.null(shiny::getDefaultReactiveDomain())) {
    shiny::showNotification(
      ui = paste(
        "Warning: Demo design file not found:",
        basename(demo_metadata_file),
        "Proceeding without experimental design information."
      ),
      id = "demo_design_file_not_found",
      duration = 10,
      type = "warning"
    )
  } else {
    warning("Demo design file not found: ", demo_metadata_file)
  }
  # Continue without design file rather than failing completely
  sample_info_demo <- NULL
} else {
  # Load design file...
}
```

### 3. General File Reading Error Handling (Lines 172-209)
**Location:** Expression file reading section

**What it does:**
- Wraps all file reading operations in `tryCatch` blocks
- Catches any file reading errors (corrupt files, permission issues, etc.)
- Provides helpful error messages to users
- Returns `NULL` gracefully on errors

**Code added:**
```r
# Wrap file reading in tryCatch for better error handling
data <- tryCatch(
  {
    # File reading logic...
  },
  error = function(e) {
    if (requireNamespace("shiny", quietly = TRUE) && !is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(
        ui = paste(
          "Error reading expression file:",
          conditionMessage(e),
          "Please check the file format and try again."
        ),
        id = "expression_file_read_error",
        duration = NULL,
        type = "error"
      )
    } else {
      warning("Error reading expression file: ", conditionMessage(e))
    }
    return(NULL)
  }
)
```

## Benefits

### User Experience
- **Clear error messages** instead of cryptic R errors
- **Helpful guidance** on what to do next
- **Graceful degradation** for missing design files (analysis can still proceed)
- **Persistent error notifications** that stay visible until dismissed

### Robustness
- **No app crashes** from missing files
- **No app crashes** from corrupt or malformed files
- **Dual context support**: Works in both Shiny and non-Shiny (testing) contexts
- **Early returns** prevent downstream errors

### Maintainability
- **Centralized error handling** in the data loading function
- **Consistent error message format** across all file operations
- **Clear comments** explaining each validation step

## Testing

A test file has been created at `tests/test_demo_file_validation.R` to verify:
1. Missing demo expression files are caught and reported
2. Missing demo design files trigger warnings but allow processing
3. File read errors are caught and reported with helpful messages

## Error Message Examples

### Missing Demo Expression File
```
Error: Demo expression file not found: missing_demo.csv
Please contact the administrator or try a different demo dataset.
```

### Missing Demo Design File
```
Warning: Demo design file not found: missing_design.csv
Proceeding without experimental design information.
```

### Corrupt/Unreadable File
```
Error reading expression file: [specific error details]
Please check the file format and try again.
```

## Files Modified
- `R/fct_01_load_data.R` - Added file validation and error handling

## Files Created
- `tests/test_demo_file_validation.R` - Test suite for new validation logic
- `DEMO_FILE_VALIDATION.md` - This documentation

## Backward Compatibility
All changes are **fully backward compatible**:
- No changes to function signatures
- No changes to expected return values for valid files
- Only adds safety checks for error conditions
