# Example script: Run NDX_Process_Subject with and without an initial spike mask
# on a synthetic dataset containing known spikes.

library(ndx)

# Generate a small synthetic dataset with injected spikes
# Similar structure to helper functions used in tests

generate_spike_dataset <- function(T_run = 30, V = 5, N_runs = 2, spike_amplitude = 50) {
  set.seed(123)
  total_T <- T_run * N_runs
  Y <- matrix(rnorm(total_T * V), nrow = total_T, ncol = V) * 0.01
  run_idx <- rep(seq_len(N_runs), each = T_run)

  # Add simple low-rank patterns
  true_V_pattern <- matrix(rnorm(V * 2), V, 2)
  for (r in seq_len(N_runs)) {
    rows <- which(run_idx == r)
    C1 <- sin((1:T_run)/8 + r/2) * 5
    C2 <- cos((1:T_run)/15 - r/3) * 4
    Y[rows, ] <- Y[rows, ] + (cbind(C1, C2) %*% t(true_V_pattern)) * 5
  }

  # Inject large spikes at known TRs
  spike_tr_list <- list(c(10, 20), c(15))
  expected_spikes_global <- integer()
  for (r in seq_len(N_runs)) {
    rows <- which(run_idx == r)
    for (tr in spike_tr_list[[r]]) {
      Y[rows[tr], ] <- spike_amplitude
      expected_spikes_global <- c(expected_spikes_global, rows[tr])
    }
  }

  list(Y_residuals_cat = Y, run_idx = run_idx,
       total_T = total_T, expected_spikes_global = expected_spikes_global)
}

# Create dataset and minimal events
TR <- 2
spike_data <- generate_spike_dataset()

motion_params <- matrix(0, nrow = spike_data$total_T, ncol = 3)
colnames(motion_params) <- paste0("mot", 1:3)

events <- data.frame(
  onsets = c(5, 15, 5, 15) * TR,
  durations = rep(2 * TR, 4),
  condition = factor(c("TaskA", "TaskB", "TaskA", "TaskB")),
  blockids = c(1, 1, 2, 2)
)

# Workflow options (reused from tests)
user_opts <- list(
  opts_pass0 = list(poly_degree = 1),
  opts_hrf = list(
    hrf_fir_taps = 6,
    hrf_fir_span_seconds = 12,
    good_voxel_R2_threshold = -Inf,
    lambda1_grid = c(0.1),
    lambda2_grid = c(0.1),
    cv_folds = 2,
    hrf_min_good_voxels = 1,
    hrf_cluster_method = "none",
    num_hrf_clusters = 1
  ),
  opts_rpca = list(k_global_target = 2, rpca_lambda_auto = FALSE, rpca_lambda_fixed = 0.1),
  opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
  opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
  opts_ridge = list(lambda_ridge = 0.5),
  max_passes = 1,
  min_des_gain_convergence = -Inf,
  min_rho_noise_projection_convergence = -Inf
)

# --- Run without an initial spike mask ---
res_no_mask <- NDX_Process_Subject(
  Y_fmri = spike_data$Y_residuals_cat,
  events = events,
  motion_params = motion_params,
  run_idx = spike_data$run_idx,
  TR = TR,
  user_options = user_opts,
  verbose = TRUE
)

# --- Run with an initial spike mask that marks the first TR ---
initial_mask <- rep(FALSE, spike_data$total_T)
initial_mask[1] <- TRUE

res_with_mask <- NDX_Process_Subject(
  Y_fmri = spike_data$Y_residuals_cat,
  events = events,
  motion_params = motion_params,
  run_idx = spike_data$run_idx,
  TR = TR,
  spike_TR_mask = initial_mask,
  user_options = user_opts,
  verbose = TRUE
)

# Compare number of valid TRs used for HRF estimation via prepare_hrf_response_data
prep_no_mask <- ndx:::prepare_hrf_response_data(
  Y_fmri = spike_data$Y_residuals_cat,
  pass0_residuals = res_no_mask$pass0_residuals,
  run_idx = spike_data$run_idx,
  validated_spike_TR_mask = res_no_mask$spike_TR_mask,
  user_options = user_opts$opts_hrf
)

prep_with_mask <- ndx:::prepare_hrf_response_data(
  Y_fmri = spike_data$Y_residuals_cat,
  pass0_residuals = res_with_mask$pass0_residuals,
  run_idx = spike_data$run_idx,
  validated_spike_TR_mask = res_with_mask$spike_TR_mask,
  user_options = user_opts$opts_hrf
)

cat(sprintf("ybar_clean length (no mask): %d\n", length(prep_no_mask$ybar_clean)))
cat(sprintf("ybar_clean length (with mask): %d\n", length(prep_with_mask$ybar_clean)))

# Confirm that the final spike mask contains the initial spikes and the detected ones
cat("Final spike indices with mask:", which(res_with_mask$spike_TR_mask), "\n")
cat("Injected spikes:", paste(spike_data$expected_spikes_global, collapse = ","), "\n")
stopifnot(all(spike_data$expected_spikes_global %in% which(res_with_mask$spike_TR_mask)))
stopifnot(all(which(initial_mask) %in% which(res_with_mask$spike_TR_mask)))

