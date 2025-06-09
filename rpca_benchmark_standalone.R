# Self-contained RPCA Benchmark
# Comparing: rpca, ADMM::admm.rpca, rsvd::rrpca

library(rpca)
library(ADMM) 
library(rsvd)

cat("=== RPCA Implementation Benchmark ===\n\n")

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

# Wrapper functions for consistent interface
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

run_rpca_rsvd <- function(M, lambda = NULL) {
  if (is.null(lambda)) lambda <- 1/sqrt(max(nrow(M), ncol(M)))
  result <- rsvd::rrpca(M, lambda = lambda)
  list(L = result$L, S = result$S)
}

# Test scenarios
scenarios <- list(
  list(name = "small", m = 40, n = 30, rank = 3, sparsity = 0.1),
  list(name = "medium", m = 60, n = 50, rank = 4, sparsity = 0.15),
  list(name = "sparse", m = 50, n = 40, rank = 3, sparsity = 0.05),
  list(name = "low_rank", m = 50, n = 40, rank = 2, sparsity = 0.1),
  list(name = "rectangular", m = 80, n = 40, rank = 4, sparsity = 0.1)
)

all_results <- data.frame()

for (scenario in scenarios) {
  cat(sprintf("Testing %s (%dx%d, rank=%d, sparsity=%.2f)...\n", 
             scenario$name, scenario$m, scenario$n, scenario$rank, scenario$sparsity))
  
  # Generate test data
  test_data <- generate_rpca_test_data(
    m = scenario$m, n = scenario$n, 
    rank = scenario$rank, sparsity_level = scenario$sparsity
  )
  
  # Test each method
  methods <- list(rpca = run_rpca_current, admm = run_rpca_admm, rsvd = run_rpca_rsvd)
  
  for (method_name in names(methods)) {
    cat(sprintf("  %s... ", method_name))
    
    # Time and test
    start_time <- Sys.time()
    tryCatch({
      result <- methods[[method_name]](test_data$M)
      end_time <- Sys.time()
      
      # Compute key metrics
      recon_error <- norm(test_data$M - (result$L + result$S), "F") / norm(test_data$M, "F")
      L_cor <- cor(c(test_data$L_true), c(result$L))
      S_cor <- cor(c(test_data$S_true), c(result$S))
      time_sec <- as.numeric(end_time - start_time)
      
      # Nuclear norm and sparsity checks
      nuclear_norm_L <- sum(svd(result$L, nu = 0, nv = 0)$d)
      nuclear_norm_true <- sum(svd(test_data$L_true, nu = 0, nv = 0)$d)
      nuclear_ratio <- nuclear_norm_L / nuclear_norm_true
      
      sparsity_recovered <- sum(abs(result$S) < 1e-6) / length(result$S)
      
      # Store result
      result_row <- data.frame(
        scenario = scenario$name,
        method = method_name,
        m = scenario$m, n = scenario$n, rank = scenario$rank, sparsity = scenario$sparsity,
        time_sec = time_sec,
        reconstruction_error = recon_error,
        L_correlation = L_cor,
        S_correlation = S_cor,
        nuclear_ratio = nuclear_ratio,
        sparsity_recovered = sparsity_recovered,
        success = TRUE,
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(all_results, result_row)
      cat(sprintf("✅ (%.2fs, recon_err=%.4f, L_cor=%.3f, S_cor=%.3f)\n", 
                 time_sec, recon_error, L_cor, S_cor))
      
    }, error = function(e) {
      # Store failed result
      result_row <- data.frame(
        scenario = scenario$name,
        method = method_name,
        m = scenario$m, n = scenario$n, rank = scenario$rank, sparsity = scenario$sparsity,
        time_sec = NA,
        reconstruction_error = Inf,
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

# Generate comprehensive summary
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("COMPREHENSIVE RESULTS SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

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

# Performance for successful runs only
successful_results <- all_results[all_results$success, ]

if (nrow(successful_results) > 0) {
  cat("\nPERFORMANCE METRICS (successful runs only):\n")
  
  # Key metrics by method
  metrics <- c("time_sec", "reconstruction_error", "L_correlation", "S_correlation", "nuclear_ratio")
  
  for (metric in metrics) {
    cat(sprintf("\n%s:\n", toupper(gsub("_", " ", metric))))
    
    metric_summary <- aggregate(successful_results[[metric]], 
                               by = list(method = successful_results$method), 
                               FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                  median = median(x, na.rm = TRUE),
                                                  min = min(x, na.rm = TRUE),
                                                  max = max(x, na.rm = TRUE)))
    
    for (i in 1:nrow(metric_summary)) {
      stats <- metric_summary$x[i, ]
      cat(sprintf("  %s: mean=%.4f, median=%.4f, range=[%.4f, %.4f]\n", 
                 metric_summary$method[i], stats["mean"], stats["median"], 
                 stats["min"], stats["max"]))
    }
  }
  
  # Overall ranking with weighted composite score
  cat("\nOVERALL RANKING:\n")
  cat("(Scoring: Quality=70%, Speed=30%)\n")
  
  # Composite score emphasizing quality over speed
  successful_results$composite_score <- (
    # Quality components (70% weight)
    -log10(successful_results$reconstruction_error + 1e-12) * 3 +     # Primary metric
    successful_results$L_correlation * 2 +                            # L component quality  
    successful_results$S_correlation * 2 +                            # S component quality
    -log10(abs(successful_results$nuclear_ratio - 1) + 1e-6) * 1 +    # Nuclear norm accuracy
    
    # Speed component (30% weight)  
    -log10(successful_results$time_sec + 0.001) * 1.5                 # Speed bonus
  )
  
  ranking <- aggregate(composite_score ~ method, successful_results, mean)
  ranking <- ranking[order(-ranking$composite_score), ]
  
  for (i in 1:nrow(ranking)) {
    cat(sprintf("  %d. %s (composite score: %.2f)\n", 
               i, ranking$method[i], ranking$composite_score[i]))
  }
  
  # Specific recommendations
  cat("\nRECOMMENDATIONS:\n")
  
  # Overall winner
  best_method <- ranking$method[1]
  cat(sprintf("🏆 OVERALL WINNER: %s\n", toupper(best_method)))
  
  # Speed ranking
  time_ranking <- aggregate(time_sec ~ method, successful_results, median)
  time_ranking <- time_ranking[order(time_ranking$time_sec), ]
  
  # Quality ranking  
  quality_ranking <- aggregate(reconstruction_error ~ method, successful_results, median)
  quality_ranking <- quality_ranking[order(quality_ranking$reconstruction_error), ]
  
  cat(sprintf("⚡ FASTEST: %s (%.3fs median time)\n", 
             time_ranking$method[1], time_ranking$time_sec[1]))
  cat(sprintf("🎯 HIGHEST QUALITY: %s (%.6f median error)\n", 
             quality_ranking$method[1], quality_ranking$reconstruction_error[1]))
  
  # Stability assessment
  error_stability <- aggregate(reconstruction_error ~ method, successful_results, sd)
  error_stability <- error_stability[order(error_stability$reconstruction_error), ]
  cat(sprintf("🔒 MOST STABLE: %s (%.6f error std dev)\n", 
             error_stability$method[1], error_stability$reconstruction_error[1]))
  
  # Use case recommendations
  cat("\nUSE CASE RECOMMENDATIONS:\n")
  cat(sprintf("• For production workflows (balance of quality & speed): %s\n", best_method))
  cat(sprintf("• For real-time applications (speed critical): %s\n", time_ranking$method[1]))
  cat(sprintf("• For research/analysis (quality critical): %s\n", quality_ranking$method[1]))
  
} else {
  cat("\n❌ No successful runs to analyze.\n")
}

# Save results
write.csv(all_results, "rpca_benchmark_results.csv", row.names = FALSE)
cat(sprintf("\n📁 Results saved to: rpca_benchmark_results.csv\n"))

# Show detailed results table
cat("\n", paste(rep("=", 100), collapse=""), "\n")
cat("DETAILED RESULTS TABLE\n")
cat(paste(rep("=", 100), collapse=""), "\n")

# Print formatted table
result_table <- all_results[, c("scenario", "method", "success", "time_sec", 
                               "reconstruction_error", "L_correlation", "S_correlation")]
result_table$time_sec <- round(result_table$time_sec, 3)
result_table$reconstruction_error <- format(result_table$reconstruction_error, scientific = TRUE, digits = 3)
result_table$L_correlation <- round(result_table$L_correlation, 3)
result_table$S_correlation <- round(result_table$S_correlation, 3)

print(result_table, row.names = FALSE)

cat("\n=== Benchmark Complete! ===\n")

# Return results for further analysis
invisible(all_results) 