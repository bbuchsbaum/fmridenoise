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