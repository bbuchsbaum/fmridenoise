# Helper to generate some Y_residuals_cat and run_idx for RPCA tests
.generate_rpca_test_data <- function(T_run = 100, V = 50, N_runs = 1) {
  total_T <- T_run * N_runs
  Y_res_cat <- matrix(rnorm(total_T * V), total_T, V) * 0.01
  run_idx_vec <- rep(1:N_runs, each = T_run)
  
  if (V >= 2) { 
      true_V_pattern <- matrix(rnorm(V*2), V, 2) 
      for (r in 1:N_runs) {
        run_rows <- which(run_idx_vec == r)
        # Stronger low-rank signal
        C_r_signal_run1 <- sin((1:T_run)/8 + r/2) * 5 
        C_r_signal_run2 <- cos((1:T_run)/15 - r/3) * 4 
        Y_res_cat[run_rows, ] <- Y_res_cat[run_rows, ] + 
                                   (cbind(C_r_signal_run1, C_r_signal_run2) %*% t(true_V_pattern)) * 5.0 # Increased signal multiplier
      }
  }
  return(list(Y_residuals_cat = Y_res_cat, run_idx = run_idx_vec, T_run=T_run, V=V, N_runs=N_runs, total_T=total_T))
} 