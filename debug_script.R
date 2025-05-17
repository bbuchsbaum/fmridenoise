# devtools::load_all("/Users/bbuchsbaum/code/fmridenoise") # Load the package

# # --- Basic Test Data Setup (copied from test-design_matrix.R) ---
# TR_test_dm <- 2.0
# n_time_per_run_dm <- 30 
# n_runs_dm <- 2
# total_timepoints_dm <- n_time_per_run_dm * n_runs_dm
# run_idx_dm <- rep(1:n_runs_dm, each = n_time_per_run_dm)

# create_mock_hrfs <- function(conditions = c("TaskA", "TaskB"), taps_per_hrf = 6, TR = 2.0) {
#   if (length(conditions) == 0) return(NULL)
#   hrf_list <- lapply(conditions, function(cond) {
#     list(
#       condition = cond,
#       hrf_estimate = list(stats::rnorm(taps_per_hrf)),
#       taps = list(1:taps_per_hrf)
#     )
#   })
#   df_hrfs <- do.call(rbind, lapply(hrf_list, function(row_list) data.frame(condition=row_list$condition)))
#   df_hrfs$hrf_estimate <- lapply(hrf_list, function(row_list) row_list$hrf_estimate[[1]])
#   df_hrfs$taps <- lapply(hrf_list, function(row_list) row_list$taps[[1]])
#   return(tibble::as_tibble(df_hrfs))
# }

# events_dm <- data.frame(
#   onsets = c(5, 15, 5, 20) * TR_test_dm,
#   durations = rep(2 * TR_test_dm, 4),
#   condition = factor(c("TaskA", "TaskB", "TaskA", "TaskB")),
#   blockids = c(1, 1, 2, 2)
# )

# motion_params_dm <- matrix(stats::rnorm(total_timepoints_dm * 3), ncol = 3)
# colnames(motion_params_dm) <- paste0("mot", 1:3)

# rpca_comps_dm <- matrix(stats::rnorm(total_timepoints_dm * 2), ncol = 2)

# spectral_sines_dm <- matrix(stats::rnorm(total_timepoints_dm * 4), ncol = 4)
# colnames(spectral_sines_dm) <- paste0(rep(c("s1", "s2"), each=2), c("_sin", "_cos"))

# estimated_hrfs_dm <- create_mock_hrfs(conditions = c("TaskA", "TaskB"), taps_per_hrf = 8, TR = TR_test_dm)
# # --- End Test Data Setup ---

# # Call the function from the first test
# message("--- Running ndx_build_design_matrix for first test case ---")
# X_full <- ndx::ndx_build_design_matrix(
#   estimated_hrfs = estimated_hrfs_dm,
#   events = events_dm,
#   motion_params = motion_params_dm,
#   rpca_components = rpca_comps_dm,
#   spectral_sines = spectral_sines_dm,
#   run_idx = run_idx_dm,
#   TR = TR_test_dm,
#   poly_degree = 1,
#   verbose = TRUE
# )

# message("--- Result for X_full ---")
# if (is.null(X_full)) {
#   message("X_full is NULL")
# } else {
#   message(sprintf("dim(X_full): %s", paste(dim(X_full), collapse="x")))
#   message(sprintf("colnames(X_full): %s", paste(colnames(X_full), collapse=", ")))
# } 


# --- Debugging for "handles single run correctly" test ---
devtools::load_all("/Users/bbuchsbaum/code/fmridenoise")

TR_test_dm <- 2.0
n_time_per_run_dm <- 30 
motion_params_dm <- matrix(stats::rnorm((n_time_per_run_dm*2) * 3), ncol = 3) # Full data size
colnames(motion_params_dm) <- paste0("mot", 1:3)
rpca_comps_dm <- matrix(stats::rnorm((n_time_per_run_dm*2) * 2), ncol = 2)    # Full data size
spectral_sines_dm <- matrix(stats::rnorm((n_time_per_run_dm*2) * 4), ncol = 4) # Full data size
colnames(spectral_sines_dm) <- paste0(rep(c("s1", "s2"), each=2), c("_sin", "_cos"))
events_dm <- data.frame( # Full events_dm needed for create_mock_hrfs if it uses it implicitly, and for subsetting
  onsets = c(5, 15, 5, 20) * TR_test_dm,
  durations = rep(2 * TR_test_dm, 4),
  condition = factor(c("TaskA", "TaskB", "TaskA", "TaskB")),
  blockids = c(1, 1, 2, 2) 
)
create_mock_hrfs <- function(conditions = c("TaskA", "TaskB"), taps_per_hrf = 6, TR = 2.0) {
  if (length(conditions) == 0) return(NULL)
  hrf_list <- lapply(conditions, function(cond) {
    list(condition = cond, hrf_estimate = list(stats::rnorm(taps_per_hrf)), taps = list(1:taps_per_hrf))
  })
  df_hrfs <- do.call(rbind, lapply(hrf_list, function(row_list) data.frame(condition=row_list$condition)))
  df_hrfs$hrf_estimate <- lapply(hrf_list, function(row_list) row_list$hrf_estimate[[1]])
  df_hrfs$taps <- lapply(hrf_list, function(row_list) row_list$taps[[1]])
  return(tibble::as_tibble(df_hrfs))
}
estimated_hrfs_dm <- create_mock_hrfs(conditions = c("TaskA", "TaskB"), taps_per_hrf = 8, TR = TR_test_dm)

# Inputs for the single run test
run_idx_single_dm <- rep(1, n_time_per_run_dm)
events_single_dm <- events_dm[events_dm$blockids == 1,]
motion_single_dm <- motion_params_dm[1:n_time_per_run_dm, , drop=FALSE]
rpca_single_dm <- rpca_comps_dm[1:n_time_per_run_dm, , drop=FALSE]
spectral_single_dm <- spectral_sines_dm[1:n_time_per_run_dm, , drop=FALSE]

message("--- Running ndx_build_design_matrix for SINGLE RUN case ---")
X_single_run <- ndx::ndx_build_design_matrix(
  estimated_hrfs = estimated_hrfs_dm, 
  events = events_single_dm,
  motion_params = motion_single_dm,
  rpca_components = rpca_single_dm,
  spectral_sines = spectral_single_dm,
  run_idx = run_idx_single_dm,
  TR = TR_test_dm,
  poly_degree = 1, 
  verbose = TRUE
)

message("--- Result for X_single_run ---")
if (is.null(X_single_run)) {
  message("X_single_run is NULL")
} else {
  message(sprintf("dim(X_single_run): %s", paste(dim(X_single_run), collapse="x")))
  message(sprintf("colnames(X_single_run): %s", paste(colnames(X_single_run), collapse=", ")))
} 

# --- Test Data Setup from test-workflow.R ---
TR_test <- 2.0
n_time_per_run_test <- 50 
n_runs_test <- 1
total_timepoints_test <- n_time_per_run_test * n_runs_test
n_voxels_test <- 5 

Y_fmri_test <- matrix(rnorm(total_timepoints_test * n_voxels_test), 
                      nrow = total_timepoints_test, 
                      ncol = n_voxels_test)
run_idx_test <- rep(1:n_runs_test, each = n_time_per_run_test)
motion_params_test <- matrix(rnorm(total_timepoints_test * 3), 
                             nrow = total_timepoints_test, ncol = 3)
colnames(motion_params_test) <- paste0("mot", 1:3)

events_test <- data.frame(
  onsets = as.numeric(c(10, 30) * TR_test), 
  durations = as.numeric(c(5, 5) * TR_test),   
  condition = factor(c("TaskA", "TaskB")),
  blockids = as.integer(rep(1, 2)) 
)
if (n_runs_test > 1) {
  events_run2_test <- data.frame(
    onsets = as.numeric(c(10, 30) * TR_test),
    durations = as.numeric(c(5, 5) * TR_test),
    condition = factor(c("TaskA", "TaskB")),
    blockids = as.integer(rep(2, 2))
  )
  events_test <- rbind(events_test, events_run2_test)
}

# User options from test-workflow.R, ensuring all HRF options are present
default_opts_hrf_from_workflow_code <- list(
    hrf_fir_taps = 12L,
    hrf_fir_span_seconds = 24, # Default from NDX_Process_Subject
    good_voxel_R2_threshold = 0.05,
    cv_folds = 5L,
    lambda1_grid = 10^seq(-2, 1, length.out = 5),
    lambda2_grid = 10^seq(-3, 0, length.out = 5),
    hrf_min_good_voxels = 50L, # Default from NDX_Process_Subject
    return_full_model = FALSE,
    hrf_cluster_method = "none", 
    num_hrf_clusters = 1 
)

user_options_test_from_workflow_file <- list(
  opts_pass0 = list(
    poly_degree = 1 
  ),
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
  opts_rpca = list(
    k_global_target = 2, 
    rpca_lambda_auto = FALSE,
    rpca_lambda_fixed = 0.1 
  ),
  opts_spectral = list(
    n_sine_candidates = 2, 
    nyquist_guard_factor = 0.1
  ),
  opts_whitening = list(
    global_ar_on_design = FALSE,
    max_ar_failures_prop = 0.5
  ),
  opts_ridge = list(
    lambda_ridge = 0.5
  ),
  task_regressor_names_for_extraction = c("task_TaskA", "task_TaskB"),
  max_passes = 2, 
  min_des_gain_convergence = -Inf, 
  min_rho_noise_projection_convergence = -Inf 
)

# Ensure all required hrf options are present, using defaults from NDX_Process_Subject if necessary
final_opts_hrf_for_debug <- utils::modifyList(default_opts_hrf_from_workflow_code, 
                                          user_options_test_from_workflow_file$opts_hrf)
user_options_for_debug_run <- user_options_test_from_workflow_file
user_options_for_debug_run$opts_hrf <- final_opts_hrf_for_debug

# --- Running NDX_Process_Subject for workflow test case ---
message("--- Running NDX_Process_Subject for WORKFLOW test case ---")

workflow_output <- ndx::NDX_Process_Subject(
  Y_fmri = Y_fmri_test,
  events = events_test,
  motion_params = motion_params_test,
  run_idx = run_idx_test,
  TR = TR_test,
  user_options = user_options_for_debug_run,
  verbose = TRUE # CRITICAL: Ensure this verbose is TRUE
)

message("--- Result for workflow_output ---")
if (is.null(workflow_output)) {
  message("workflow_output is NULL")
} else {
  message(sprintf("workflow_output names: %s", paste(names(workflow_output), collapse=", ")))
  message(sprintf("Number of passes completed: %d", workflow_output$num_passes_completed))
  if (length(workflow_output$diagnostics_per_pass) > 0) {
      message("Diagnostics for pass 1 DES:")
      print(workflow_output$diagnostics_per_pass[[1]]$DES)
      message("Diagnostics for pass 1 Rho:")
      print(workflow_output$diagnostics_per_pass[[1]]$rho_noise_projection)
  }
} 