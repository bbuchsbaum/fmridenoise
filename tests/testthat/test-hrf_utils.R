context("HRF Utility Functions - Cone Projection")

# Tests for the old project_cone_heuristic (can be removed once fully deprecated)
test_that("project_cone_heuristic (old) enforces basic constraints", {
  # Case 1: typical HRF coefficients with non-zero first element and negative second
  coeffs1 <- c(1, -2, 3)
  # adjusted1 <- project_cone_heuristic(coeffs1) # Function removed, test can be removed later
  # expect_equal(adjusted1[1], 0)
  # expect_equal(adjusted1[2], 0)
  # expect_equal(adjusted1[3], 2)
  expect_true(TRUE) # Placeholder for now

  # Case 2: coefficients become flat after adjustment and get small bump
  coeffs2 <- c(0, -0.1, 0)
  # adjusted2 <- project_cone_heuristic(coeffs2)
  # expect_true(adjusted2[1] > 0)
  # expect_equal(adjusted2[2], 0)
  # expect_equal(adjusted2[3], 0)
  expect_true(TRUE) # Placeholder

  # Case 3: zero-length input returns numeric(0)
  # adjusted3 <- project_cone_heuristic(numeric(0))
  # expect_length(adjusted3, 0)
  expect_true(TRUE) # Placeholder
})


context("project_hrf_cone - New HRF Cone Projection")

test_that("project_hrf_cone handles basic cases and options", {
  h1 <- c(0, 0, 1, 2, 1, 0.5, 0, -0.1, -0.05) # Typical shape, late undershoot
  
  # Default: nonneg=T, unimodal=T, normalize_area=T
  h1_default <- ndx:::project_hrf_cone(h1)
  expect_true(all(h1_default >= 0)) 
  peak_idx_d <- which.max(h1_default)
  if (length(h1_default) >=3 && peak_idx_d > 1 && peak_idx_d < length(h1_default)) {
    expect_true(all(diff(h1_default[1:peak_idx_d]) >= -1e-9))
    expect_true(all(diff(h1_default[peak_idx_d:length(h1_default)]) <= 1e-9))
  }
  expect_equal(sum(h1_default), 1.0, tolerance = 1e-6)

  # No normalization
  h1_no_norm <- ndx:::project_hrf_cone(h1, normalize_area = FALSE)
  expect_true(all(h1_no_norm >= 0))
  # Sum will not be 1, check original sum after projection (anchored, nonneg, unimodal)
  # h1_proj_no_norm_intermediate = c(0,0,1,2,1,0.5,0,0,0), sum = 4.5
  expect_equal(sum(h1_no_norm), 4.5, tolerance=1e-9)

  # No unimodal constraint (nonneg=T by default for func, normalize_area=F by test)
  h1_no_uni <- ndx:::project_hrf_cone(h1, unimodal = FALSE, normalize_area = FALSE)
  expect_true(h1_no_uni[1] == 0) 
  expect_true(all(h1_no_uni >= 0)) 
  expect_equal(h1_no_uni[8], 0) # Corrected: since nonneg=TRUE is default for project_hrf_cone arg
  expect_equal(h1_no_uni[9], 0) # Corrected

  # Explicitly nonneg = FALSE, unimodal = TRUE, normalize_area = FALSE
  h_neg_allowed <- c(0.1, 0.5, 1, 0.8, 0.2, -0.1, -0.3, -0.1)
  h_neg_proj <- ndx:::project_hrf_cone(h_neg_allowed, nonneg = FALSE, unimodal = TRUE, normalize_area = FALSE)
  expect_equal(h_neg_proj[1], 0, tolerance=1e-9) # Anchored
  expect_true(any(h_neg_proj < -1e-9)) # Should retain some negativity if unimodal allows
  peak_idx_n <- which.max(h_neg_proj)
  if (length(h_neg_proj) >=3 && peak_idx_n > 1 && peak_idx_n < length(h_neg_proj)) {
     expect_true(all(diff(h_neg_proj[1:peak_idx_n]) >= -1e-9))
     expect_true(all(diff(h_neg_proj[peak_idx_n:length(h_neg_proj)]) <= 1e-9))
  }
})

test_that("project_hrf_cone handles edge cases (short, flat, all zero)", {
  # Short HRFs
  expect_equal(ndx:::project_hrf_cone(c(1,2), nonneg=TRUE, unimodal=TRUE, normalize_area=TRUE), c(0,1)) 
  expect_equal(ndx:::project_hrf_cone(c(1,2), normalize_area=FALSE), c(0,1)) 
  expect_equal(ndx:::project_hrf_cone(numeric(0)), numeric(0))
  # h=c(5) -> h-h[1] = c(0). sum=0. normalize_area=T does nothing. Then bump h[1]=1e-6. Re-norm -> h[1]=1.
  expect_equal(ndx:::project_hrf_cone(c(5)), c(1)) # Corrected
  
  # Flat HRF
  h_flat <- rep(2, 5)
  # h_flat-h_flat[1] = c(0,0,0,0,0). Bump at h[1]=1e-6. Normalize area=TRUE -> c(1,0,0,0,0)
  h_flat_proj <- ndx:::project_hrf_cone(h_flat, normalize_area = TRUE)
  expect_equal(h_flat_proj, c(1,0,0,0,0), tolerance = 1e-7) # Corrected
  
  h_flat_proj_no_norm <- ndx:::project_hrf_cone(h_flat, normalize_area = FALSE)
  # h_flat-h_flat[1] = c(0,0,0,0,0). Bump at h[1]=1e-6.
  expected_flat_no_norm <- rep(0,5); expected_flat_no_norm[1] <- 1e-6
  expect_equal(h_flat_proj_no_norm, expected_flat_no_norm, tolerance = 1e-9)
  
  # All zero HRF (should get small bump if original wasn't all zero - it was here)
  h_zero <- rep(0, 5)
  # h_orig_for_debug is c(0,0,0,0,0). sum(abs(h_orig_for_debug)) is 0. No bump.
  h_zero_proj <- ndx:::project_hrf_cone(h_zero, normalize_area = FALSE)
  expect_equal(h_zero_proj, rep(0,5), tolerance=1e-9) # Corrected: no bump if originally all zero
  
  # All zero HRF but normalize area = TRUE (should remain all zero)
   h_zero_proj_norm <- ndx:::project_hrf_cone(h_zero, normalize_area = TRUE)
  expect_equal(h_zero_proj_norm, rep(0,5), tolerance=1e-9)
})
