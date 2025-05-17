context("spike mask validation")

# Use internal helper from ndx_rpca.R


test_that("finalize_spike_mask warns and fixes length mismatch", {
  mask_list <- list(run_1 = c(TRUE, FALSE),
                    run_2 = c(FALSE))
  unique_runs <- c(1, 2)
  n_total <- 5
  expect_warning(
    mask_final <- .finalize_spike_mask(mask_list, unique_runs, n_total),
    "spike_TR_mask length"
  )
  expect_length(mask_final, n_total)
  expect_type(mask_final, "logical")
})
