# Final RPCA Implementation Comparison
# Comparing: ADMM, RPCA (current), and Optimized RSVD (rsvd_premium)

library(rpca)
library(ADMM) 
library(rsvd)

cat("=== FINAL RPCA IMPLEMENTATION COMPARISON ===\n")
cat("Optimized RSVD vs. Baseline Methods\n\n")

# Generate synthetic test data with known ground truth
generate_rpca_test_data <- function(m = 100, n = 80, rank = 5, 
                                   sparsity_level = 0.1, noise_level = 0.01,
                                   sparse_magnitude = 1.0) {
  
  # Generate low-rank component L = UV^T
  U <- matrix(rnorm(m * rank), m, rank)
  V <- matrix(rnorm(n * rank), n, rank)
  L_true <- U %*% t(V)
  
  # Generate sparse component S
  S_true <- matrix(0, m, n)
  n_sparse <- round(sparsity_level * m * n)
  sparse_indices <- sample(m * n, n_sparse)
  S_true[sparse_indices] <- rnorm(n_sparse, 0, sparse_magnitude)
  
  # Add noise
  noise <- matrix(rnorm(m * n, 0, noise_level), m, n)
  
  # Observed matrix
  M <- L_true + S_true + noise
  
  list(M = M, L_true = L_true, S_true = S_true, noise = noise,
       rank = rank, sparsity_level = sparsity_level)
}

# Method implementations
run_rpca_current <- function(M, lambda = NULL) {
  if (is.null(lambda)) lambda <- 1/sqrt(max(nrow(M), ncol(M)))
  result <- rpca::rpca(M, lambda = lambda)
  list(L = result$L, S = result$S)
}

run_rpca_admm <- function(M, lambda = NULL) {
  if (is.null(lambda)) lambda <- 1/sqrt(max(nrow(M), ncol(M)))
  result <- ADMM::admm.rpca(M, lambda = lambda)
  list(L = result$L, S = result$S)
}

# OPTIMIZED RSVD - Winner configuration from parameter tuning!
run_rsvd_optimized <- function(M, lambda = NULL) {
  if (is.null(lambda)) lambda <- 1/sqrt(max(nrow(M), ncol(M)))
  # Optimal parameters: maxiter=300, tol=1e-09, p=25, q=5
  result <- rsvd::rrpca(M, lambda = lambda, maxiter = 300, tol = 1e-09, p = 25, q = 5)
  list(L = result$L, S = result$S)
}

# Test scenarios - comprehensive coverage
scenarios <- list(
  list(name = "small", m = 40, n = 30, rank = 3, sparsity = 0.1),
  list(name = "medium", m = 60, n = 50, rank = 4, sparsity = 0.15),
  list(name = "large", m = 100, n = 80, rank = 6, sparsity = 0.1),
  list(name = "sparse", m = 50, n = 40, rank = 3, sparsity = 0.05),
  list(name = "dense", m = 50, n = 40, rank = 3, sparsity = 0.25),
  list(name = "low_rank", m = 60, n = 50, rank = 2, sparsity = 0.1),
  list(name = "high_rank", m = 60, n = 50, rank = 10, sparsity = 0.1),
  list(name = "rectangular_wide", m = 40, n = 80, rank = 4, sparsity = 0.1),
  list(name = "rectangular_tall", m = 80, n = 40, rank = 4, sparsity = 0.1)
)

all_results <- data.frame()

for (scenario in scenarios) {
  cat(sprintf("\nTesting %s (%dx%d, rank=%d, sparsity=%.2f)...\n", 
             scenario$name, scenario$m, scenario$n, scenario$rank, scenario$sparsity))
  
  # Generate test data
  test_data <- generate_rpca_test_data(
    m = scenario$m, n = scenario$n, 
    rank = scenario$rank, sparsity_level = scenario$sparsity
  )
  
  # Test each method
  methods <- list(
    list(name = "rpca_current", func = run_rpca_current),
    list(name = "admm", func = run_rpca_admm),
    list(name = "rsvd_optimized", func = run_rsvd_optimized)
  )
  
  for (method in methods) {
    cat(sprintf("  %s... ", method$name))
    
    start_time <- Sys.time()
    tryCatch({
      result <- method$func(test_data$M)
      end_time <- Sys.time()
      
      # Compute comprehensive metrics
      recon_error <- norm(test_data$M - (result$L + result$S), "F") / norm(test_data$M, "F")
      L_cor <- cor(c(test_data$L_true), c(result$L))
      S_cor <- cor(c(test_data$S_true), c(result$S))
      time_sec <- as.numeric(end_time - start_time)
      
      # Additional quality metrics
      L_error <- norm(result$L - test_data$L_true, "F") / norm(test_data$L_true, "F")
      S_error <- norm(result$S - test_data$S_true, "F") / norm(test_data$S_true, "F")
      
      nuclear_norm_L <- sum(svd(result$L, nu = 0, nv = 0)$d)
      nuclear_norm_true <- sum(svd(test_data$L_true, nu = 0, nv = 0)$d)
      nuclear_ratio <- nuclear_norm_L / nuclear_norm_true
      
      sparsity_recovered <- sum(abs(result$S) < 1e-6) / length(result$S)
      
      # Store result
      result_row <- data.frame(
        scenario = scenario$name,
        method = method$name,
        m = scenario$m, n = scenario$n, rank = scenario$rank, sparsity = scenario$sparsity,
        time_sec = time_sec,
        reconstruction_error = recon_error,
        L_error = L_error,
        S_error = S_error,
        L_correlation = L_cor,
        S_correlation = S_cor,
        nuclear_ratio = nuclear_ratio,
        sparsity_recovered = sparsity_recovered,
        success = TRUE,
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(all_results, result_row)
      cat(sprintf("✅ (%.3fs, recon=%.2e, L_cor=%.3f, S_cor=%.3f)\n", 
                 time_sec, recon_error, L_cor, S_cor))
      
    }, error = function(e) {
      # Store failed result
      result_row <- data.frame(
        scenario = scenario$name,
        method = method$name,
        m = scenario$m, n = scenario$n, rank = scenario$rank, sparsity = scenario$sparsity,
        time_sec = NA,
        reconstruction_error = Inf,
        L_error = Inf,
        S_error = Inf,
        L_correlation = 0,
        S_correlation = 0,
        nuclear_ratio = Inf,
        sparsity_recovered = 0,
        success = FALSE,
        stringsAsFactors = FALSE
      )
      all_results <<- rbind(all_results, result_row)
      cat(sprintf("❌ FAILED: %s\n", e$message))
    })
  }
}

# Comprehensive analysis
cat("\n", paste(rep("=", 90), collapse=""), "\n")
cat("FINAL COMPARISON RESULTS - OPTIMIZED RSVD vs BASELINES\n")
cat(paste(rep("=", 90), collapse=""), "\n\n")

# Filter successful results
successful_results <- all_results[all_results$success, ]

if (nrow(successful_results) > 0) {
  
  # Success rates
  cat("SUCCESS RATES:\n")
  success_summary <- aggregate(success ~ method, all_results, function(x) sum(x) / length(x))
  for (i in 1:nrow(success_summary)) {
    cat(sprintf("  %s: %.1f%% (%d/%d scenarios)\n", 
               success_summary$method[i], 
               success_summary$success[i] * 100,
               sum(all_results$method == success_summary$method[i] & all_results$success),
               sum(all_results$method == success_summary$method[i])))
  }
  
  # Performance summary
  cat("\nPERFORMANCE SUMMARY:\n")
  
  metrics <- c("time_sec", "reconstruction_error", "L_error", "S_error", 
               "L_correlation", "S_correlation", "nuclear_ratio")
  
  perf_summary <- aggregate(successful_results[, metrics], 
                           by = list(method = successful_results$method), 
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                              median = median(x, na.rm = TRUE),
                                              min = min(x, na.rm = TRUE),
                                              max = max(x, na.rm = TRUE)))
  
  cat(sprintf("%-18s %12s %14s %10s %10s %10s %10s\n", 
             "Method", "Time (s)", "Recon Error", "L Error", "S Error", "L Corr", "S Corr"))
  cat(paste(rep("-", 90), collapse=""), "\n")
  
  for (i in 1:nrow(perf_summary)) {
    method <- perf_summary$method[i]
    cat(sprintf("%-18s %11.3f %13.2e %9.2e %9.2e %9.3f %9.3f\n", 
               toupper(method),
               perf_summary$time_sec[i, "median"],
               perf_summary$reconstruction_error[i, "median"],
               perf_summary$L_error[i, "median"],
               perf_summary$S_error[i, "median"],
               perf_summary$L_correlation[i, "median"],
               perf_summary$S_correlation[i, "median"]))
  }
  
  # Speed comparison
  cat("\nSPEED COMPARISON (vs slowest):\n")
  time_summary <- aggregate(time_sec ~ method, successful_results, median)
  time_summary$speedup <- max(time_summary$time_sec) / time_summary$time_sec
  time_summary <- time_summary[order(-time_summary$speedup), ]
  
  for (i in 1:nrow(time_summary)) {
    cat(sprintf("  %s: %.3fs (%.1fx speedup)\n", 
               toupper(time_summary$method[i]), 
               time_summary$time_sec[i], 
               time_summary$speedup[i]))
  }
  
  # Quality comparison (vs worst)
  cat("\nQUALITY COMPARISON (reconstruction error):\n")
  error_summary <- aggregate(reconstruction_error ~ method, successful_results, median)
  error_summary$improvement <- max(error_summary$reconstruction_error) / error_summary$reconstruction_error
  error_summary <- error_summary[order(-error_summary$improvement), ]
  
  for (i in 1:nrow(error_summary)) {
    cat(sprintf("  %s: %.2e error (%.1fx better than worst)\n", 
               toupper(error_summary$method[i]), 
               error_summary$reconstruction_error[i], 
               error_summary$improvement[i]))
  }
  
  # Overall ranking
  cat("\nOVERALL RANKING (weighted composite score):\n")
  cat("Weights: Quality=60%, Speed=25%, Stability=15%\n")
  
  # Calculate composite scores for each method
  method_scores <- aggregate(
    cbind(reconstruction_error, L_correlation, S_correlation, time_sec) ~ method, 
    successful_results, 
    function(x) c(
      mean = mean(x, na.rm = TRUE), 
      median = median(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  )
  
  # Composite score calculation
  composite_scores <- data.frame(method = method_scores$method)
  
  # Quality component (60% weight)
  quality_score <- (
    -log10(method_scores$reconstruction_error[, "median"] + 1e-12) * 3 +    # Primary quality
    method_scores$L_correlation[, "median"] * 1.5 +                         # L component
    method_scores$S_correlation[, "median"] * 1.5                           # S component
  )
  
  # Speed component (25% weight)  
  max_time <- max(method_scores$time_sec[, "median"])
  speed_score <- log10(max_time / method_scores$time_sec[, "median"]) * 2.5
  
  # Stability component (15% weight)
  stability_score <- -log10(method_scores$reconstruction_error[, "sd"] + 1e-12) * 1.5
  
  composite_scores$total_score <- quality_score + speed_score + stability_score
  composite_scores <- composite_scores[order(-composite_scores$total_score), ]
  
  for (i in 1:nrow(composite_scores)) {
    cat(sprintf("  %d. %s (score: %.2f)\n", 
               i, toupper(composite_scores$method[i]), composite_scores$total_score[i]))
  }
  
  # Final recommendation
  cat("\n🏆 FINAL RECOMMENDATION 🏆\n")
  winner <- composite_scores$method[1]
  cat(sprintf("WINNER: %s\n", toupper(winner)))
  
  if (winner == "rsvd_optimized") {
    cat("\n🎯 BREAKTHROUGH: Optimized RSVD is the clear winner!\n")
    cat("Parameters: maxiter=300, tol=1e-09, p=25, q=5\n")
    
    # Show improvement over current method
    current_time <- median(successful_results[successful_results$method == "rpca_current", "time_sec"])
    winner_time <- median(successful_results[successful_results$method == winner, "time_sec"])
    speedup <- current_time / winner_time
    
    current_error <- median(successful_results[successful_results$method == "rpca_current", "reconstruction_error"])
    winner_error <- median(successful_results[successful_results$method == winner, "reconstruction_error"])
    quality_ratio <- current_error / winner_error
    
    cat(sprintf("Improvement over current rpca: %.1fx faster, %.1fx better quality\n", 
               speedup, quality_ratio))
  } else {
    cat(sprintf("\nBest method: %s\n", toupper(winner)))
  }
  
} else {
  cat("No successful results to analyze.\n")
}

# Save results
write.csv(all_results, "final_rpca_comparison_results.csv", row.names = FALSE)
cat(sprintf("\n📁 Detailed results saved to: final_rpca_comparison_results.csv\n"))

cat("\n=== Final Comparison Complete! ===\n")

# Return results
invisible(all_results) 