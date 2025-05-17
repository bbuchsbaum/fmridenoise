library(testthat)
library(fmrireg) # Explicitly load fmrireg
library(stats)   # Explicitly load stats (though often default)
library(rpca) # Ensure rpca is loaded

# Get the directory of the script itself
script_dir <- getwd()

# Define paths relative to the script's location (workspace root)
R_dir <- file.path(script_dir, "R")
test_dir_path <- file.path(script_dir, "tests", "testthat")

# Source all .R files in the R directory
cat("Sourcing files from:", R_dir, "\n")
for (r_file in list.files(R_dir, pattern = "\\.R$", full.names = TRUE)) {
  cat("Sourcing:", r_file, "\n")
  tryCatch({
    source(r_file, chdir = TRUE) # chdir=TRUE so relative paths inside sourced files might work
  }, error = function(e) {
    cat("Error sourcing", r_file, ":", e$message, "\n")
    stop("Failed to source R files.") # Stop if sourcing fails
  })
}

# Check if test directory exists
if (dir.exists(test_dir_path)) {
  cat("Running tests from:", test_dir_path, "\n")
  
  # test_dir will error if any test fails when stop_on_failure = TRUE (default is FALSE for test_dir)
  # However, for scripting, it's often better to use test_file or devtools::test which handle exit codes.
  # For now, let's use tryCatch to see the behavior. A non-zero exit from Rscript would be the CI signal.
  test_outcome <- NULL
  tryCatch({
    # Use a reporter that's less verbose but clear, e.g., "check" or "minimal"
    # The default reporter for test_dir is "summary"
    # Using stop_on_failure=TRUE would make the script exit on first failure.
    # For now, let's capture the result of the default summary reporter.
    test_results_output <- test_dir(test_dir_path, reporter = "summary") 
    # If test_dir completes without error when stop_on_failure is FALSE (default),
    # we then need to parse its output or the returned object to check for failures.
    # The object returned by test_dir with reporter="summary" is a list of data frames.
    
    # Check if any test failed based on the structure of test_results_output
    # Each element in the list is a data.frame for a test file.
    # We look for any 'failed' column > 0 or 'error' column == TRUE.
    any_failures <- FALSE
    if (length(test_results_output) > 0) {
        for (res_df in test_results_output) {
            if (is.data.frame(res_df) && "failed" %in% colnames(res_df) && sum(res_df$failed) > 0) {
                any_failures <- TRUE
                break
            }
            if (is.data.frame(res_df) && "error" %in% colnames(res_df) && any(res_df$error)) {
                any_failures <- TRUE
                break
            }
        }
    }

    if (any_failures) {
        cat("\nTESTS FAILED.\n")
        q(status = 1) # Exit with status 1 for failure
    } else {
        cat("\nALL TESTS PASSED.\n")
        # q(status = 0) # Default exit is 0
    }
    
  }, error = function(e) {
    # This error is if test_dir itself fails catastrophically, not just a test expectation failure
    cat("\nError during test execution phase:", e$message, "\n")
    q(status = 2)
  })

} else {
  cat("Test directory not found:", test_dir_path, "\n")
  q(status = 3)
} 