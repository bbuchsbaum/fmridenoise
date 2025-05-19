# Helper to generate spike data for multiple runs
.generate_spike_data <- function(T_run = 30, V = 20, N_runs = 2) {
  set.seed(123)
  base <- .generate_rpca_test_data(T_run = T_run, V = V, N_runs = N_runs)
  spike_amplitude <- 50
  spike_tr_list <- list(c(10, 20), c(15))
  for (r in seq_len(N_runs)) {
    rows <- which(base$run_idx == r)
    spikes <- spike_tr_list[[r]]
    for (tr in spikes) {
      base$Y_residuals_cat[rows[tr], ] <- spike_amplitude
    }
  }
  base$expected_spikes_global <- c(10, 20, T_run + 15)
  base
}
