context("estimate_hrf_for_condition fused-lasso accuracy")

set.seed(123)

# Step 1: create a simple HRF shape
TR <- 1
hrf_length <- 8
true_hrf <- dgamma(seq(0, hrf_length-1) * TR, shape = 6, rate = 1)
true_hrf <- true_hrf / max(true_hrf)

# Step 2: events and sampling frame
n_time <- 80
run_idx <- rep(1, n_time)

events_df <- data.frame(
  onsets = c(10, 30, 50),
  durations = rep(1, 3),
  condition = "A",
  blockids = 1
)

sf <- fmrireg::sampling_frame(blocklens = n_time, TR = TR)

# Step 3: FIR design matrix and synthetic data
X <- get_fir_design_matrix_for_condition(
  condition_name = "A",
  events_df = events_df,
  sampling_frame = sf,
  fir_taps = hrf_length,
  TR = TR,
  verbose = FALSE
)

noise <- as.numeric(stats::arima.sim(list(ar = c(0.4, -0.2)), n = n_time, sd = 0.1))
ybar_clean <- as.numeric(X %*% true_hrf) + noise

# Step 4: estimate HRF using fused lasso
user_opts <- list(
  hrf_fir_taps = hrf_length,
  lambda1_grid = c(0.01, 0.1),
  lambda2_grid = c(0.001, 0.01),
  cv_folds = 2,
  hrf_min_events_for_fir = 1,
  hrf_low_event_threshold = 1,
  hrf_target_event_count_for_lambda_scaling = 20,
  hrf_use_canonical_fallback_for_ultra_sparse = FALSE,
  verbose_hrf = FALSE,
  return_full_model = FALSE,
  hrf_cone_nonneg = TRUE,
  hrf_cone_unimodal = TRUE,
  hrf_cone_normalize_area = TRUE
)

est_res <- ndx:::estimate_hrf_for_condition(
  condition_name = "A",
  events_for_condition = events_df,
  ybar_clean = ybar_clean,
  block_ids_for_cv = run_idx,
  overall_sampling_frame = sf,
  valid_TRs_mask = rep(TRUE, n_time),
  TR = TR,
  user_options = user_opts
)

# Step 5: verify accuracy and structure
true_hrf_adj <- ndx:::project_hrf_cone(true_hrf, 
                                     nonneg = user_opts$hrf_cone_nonneg %||% TRUE, 
                                     unimodal = user_opts$hrf_cone_unimodal %||% TRUE, 
                                     normalize_area = user_opts$hrf_cone_normalize_area %||% TRUE,
                                     verbose = user_opts$verbose_hrf %||% FALSE
                                     )

cor_val <- stats::cor(est_res$hrf_estimate, true_hrf_adj)


test_that("fused-lasso HRF estimate correlates with ground truth", {
  expect_true(is.list(est_res))
  expect_equal(est_res$condition, "A")
  expect_equal(length(est_res$hrf_estimate), hrf_length)
  expect_gt(cor_val, 0.8)
  expect_equal(est_res$taps, seq_len(hrf_length))
})

