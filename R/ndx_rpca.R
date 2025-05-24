#' Perform Multi-Run Robust PCA and Extract Temporal Nuisance Components
#'
#' This function implements a multi-run RPCA strategy. It performs RPCA on each
#' run's transposed residuals (voxels x time) to obtain voxel-space principal components (Vr).
#' These Vr components are then merged across runs using Grassmann averaging to find a
#' global voxel-space nuisance basis (V_global). Finally, per-run residuals are
#' projected onto V_global to get run-specific temporal nuisance regressors (Cr),
#' which are then concatenated.
#' This approach is designed to be memory-efficient and respect fMRI geometry.
#'
#' @param Y_residuals_cat A numeric matrix of concatenated residuals from all runs
#'   (total_timepoints x voxels).
#' @param run_idx A numeric vector indicating run membership for each row (timepoint)
#'   in `Y_residuals_cat`.
#' @param k_global_target Integer, the target number of global temporal nuisance
#'   components to be returned (columns in the output matrix).
#' @param user_options A list of user-configurable options:
#'   - `k_per_run_target` (integer): Target rank for the L component in per-run RPCA. 
#'     Defaults to `k_global_target` if not specified or if larger.
#'     It's the `k` passed to `rpca::rpca` for each run.
#'   - `rpca_term_delta` (numeric): Convergence tolerance for rpca (passed as `term.delta`). Default: 1e-6.
#'   - `rpca_max_iter` (integer): Maximum iterations for rpca (passed as `max.iter`). Default: 2000.
#'   - `rpca_trace` (logical): If TRUE, rpca will print progress messages. Default: FALSE.
#'   - `rpca_mu` (numeric): Mu parameter for rpca. Default: NULL (auto by rpca package).
#'   - `rpca_lambda_auto` (logical): If TRUE (default), calculate lambda for each run as 
#'     `1/sqrt(max(dim(Er_t)))`. If FALSE, use `rpca_lambda_fixed`.
#'   - `rpca_lambda_fixed` (numeric): Fixed lambda value if `rpca_lambda_auto` is FALSE.
#'   - `rpca_merge_strategy` (character): Strategy for merging voxel-space components.
#'     Options are "concat_svd" or "iterative". Default: "concat_svd".
#'   - `rpca_spike_mad_thresh` (numeric): MAD multiplier for spike detection. Default: 3.0.
#'   - `rpca_spike_percentile_thresh` (numeric): Percentile threshold if MAD is tiny. Default: 0.98.
#' @return A list with elements:
#'   - `C_components`: Matrix of concatenated temporal nuisance components (total_timepoints x k_global_target).
#'   - `spike_TR_mask`: Logical vector (length total_timepoints) flagging TRs with non-zero sparse activity.
#'   - `S_matrix_cat`: Matrix of concatenated S_r matrices (total_timepoints x voxels).
#'   - `V_global_singular_values`: Singular values from V_global if using concat_svd strategy.
#'   Returns NULL if errors occur or no components generated.
#' @examples
#' \dontrun{
#' # --- Simulate multi-run data ---
#' T_run <- 50; V <- 30; N_runs <- 2
#' total_T <- T_run * N_runs
#' Y_res_cat <- matrix(rnorm(total_T * V), total_T, V)
#' run_idx_vec <- rep(1:N_runs, each = T_run)
#' 
#' # Add some shared low-rank structure (voxel-space pattern, different temporal expression)
#' true_V_pattern <- matrix(rnorm(V*2), V, 2) # 2 global voxel patterns
#' for (r in 1:N_runs) {
#'   run_rows <- which(run_idx_vec == r)
#'   # Run-specific temporal modulation of these voxel patterns
#'   C_r_signal_run1 <- sin((1:T_run)/5 + r) * 3 + cos((1:T_run)/10 - r/2) * 2
#'   C_r_signal_run2 <- cos((1:T_run)/3 - r) * 2.5 + sin((1:T_run/8) + r/3) * 3
#'   Y_res_cat[run_rows, ] <- Y_res_cat[run_rows, ] + 
#'                            cbind(C_r_signal_run1, C_r_signal_run2) %*% t(true_V_pattern)
#' }
#' 
#' k_target_final <- 3
#' user_opts_mrpca <- list(
#'   k_per_run_target = 5, # Keep a bit more per run initially
#'   rpca_term_delta = 1e-4, # Relax tolerance for example speed
#'   rpca_max_iter = 50, # Reduced for example speed
#'   rpca_lambda_auto = TRUE
#' )
#' 
#' C_components <- ndx_rpca_temporal_components_multirun(
#'   Y_residuals_cat = Y_res_cat, 
#'   run_idx = run_idx_vec, 
#'   k_global_target = k_target_final, 
#'   user_options = user_opts_mrpca
#' )
#' 
#' if (!is.null(C_components)) {
#'   print(paste("Dimensions of concatenated C components:", 
#'               paste(dim(C_components), collapse="x")))
#'   # plot(C_components[,1], type='l', main="First Global RPCA Temporal Component")
#' }
#' }
#' @importFrom rpca rpca
#' @import stats
#' @export
ndx_rpca_temporal_components_multirun <- function(Y_residuals_cat, run_idx, 
                                                k_global_target, user_options = list()) {

  # --- 1. Input Validation & Options --- 
  if (!is.matrix(Y_residuals_cat) || !is.numeric(Y_residuals_cat)) {
    stop("Y_residuals_cat must be a numeric matrix (total_timepoints x voxels).")
  }
  if (!is.numeric(run_idx) || length(run_idx) != nrow(Y_residuals_cat)) {
    stop("run_idx length must match nrow(Y_residuals_cat).")
  }
  if (k_global_target < 0) {
    stop("k_global_target must be non-negative.")
  }
  if (k_global_target == 0) {
    message("k_global_target is 0, returning NULL as no components requested.")
    return(NULL)
  }

  default_opts <- list(
    k_per_run_target = k_global_target,
    rpca_term_delta = 1e-6,
    rpca_max_iter = 2000,
    rpca_mu = NULL,
    rpca_lambda_auto = TRUE,
    rpca_lambda_fixed = NULL,
    rpca_trace = FALSE,
    rpca_merge_strategy = "concat_svd", # "concat_svd" or "iterative"
    rpca_spike_mad_thresh = 3.0,
    rpca_spike_percentile_thresh = 0.98
  )
  current_opts <- utils::modifyList(default_opts, user_options)
  
  # Ensure k_per_run_target is reasonable
  if (current_opts$k_per_run_target < k_global_target && current_opts$k_per_run_target > 0) {
      message(sprintf("k_per_run_target (%d) is less than k_global_target (%d). Will use k_per_run_target=%d for per-run RPCA.",
                      current_opts$k_per_run_target, k_global_target, current_opts$k_per_run_target))
  } else if (current_opts$k_per_run_target <= 0) {
       message(sprintf("k_per_run_target (%d) is <=0. Setting to k_global_target (%d) for per-run RPCA.",
                      current_opts$k_per_run_target, k_global_target))
      current_opts$k_per_run_target <- k_global_target
  }
  
  # Split concatenated residuals into a list of per-run matrices
  unique_runs <- sort(unique(run_idx))
  if (length(unique_runs) == 0) {
    stop("No runs found in run_idx.")
  }
  
  Y_residuals_list <- lapply(unique_runs, function(r_id) {
    Y_residuals_cat[run_idx == r_id, , drop = FALSE]
  })
  names(Y_residuals_list) <- paste0("run_", unique_runs)

  if (length(Y_residuals_list) == 0) {
      warning("Splitting Y_residuals_cat by run_idx resulted in an empty list.")
      return(NULL)
  }
  
  # --- 2. Per-Run RPCA (on E_r^T) --- 
  V_list <- list() # To store V_r (voxel-space components from each run)
  glitch_ratios_per_run <- numeric(length(Y_residuals_list))
  names(glitch_ratios_per_run) <- names(Y_residuals_list)
  per_run_spike_TR_masks <- vector("list", length(Y_residuals_list))
  names(per_run_spike_TR_masks) <- names(Y_residuals_list)
  S_matrix_list_per_run_TpV <- list() # To store S_r (Time x Voxels) for each run
  V_global_singular_values <- NULL # To store singular values for rank adaptation

  message(sprintf("Starting per-run RPCA for %d runs...", length(Y_residuals_list)))
  for (r_idx in seq_along(Y_residuals_list)) {
    run_name <- names(Y_residuals_list)[r_idx]
    Er <- Y_residuals_list[[r_idx]] # Time_r x Voxels
    
    # Initialize spike mask for this run to all FALSE.
    # This ensures it's populated even if later steps fail for this run.
    # Also ensures it has an entry if a run is skipped early.
    if (nrow(Er) > 0) {
        per_run_spike_TR_masks[[run_name]] <- rep(FALSE, nrow(Er))
    } else {
        # If Er is empty, the mask should also be empty for unlist to work correctly later
        # Or handle this case specifically where Y_residuals_list has 0-row entries.
        # For now, if Er has 0 rows, this run will likely be skipped anyway.
        # Let's assume nrow(Er) > 0 if we reach here for spike mask generation.
        # The existing check below handles empty Er for RPCA itself.
    }
    
    if (nrow(Er) == 0 || ncol(Er) == 0) {
        warning(sprintf("Residuals for run %s are empty (dims: %s). Skipping RPCA for this run.", 
                        run_name, paste(dim(Er), collapse="x")))
        V_list[[run_name]] <- NULL # Placeholder for potential later filtering
        glitch_ratios_per_run[run_name] <- NA
        # per_run_spike_TR_masks[[run_name]] already set to FALSE vector of length nrow(Er) if nrow(Er) > 0
        # or remains NULL if nrow(Er) was 0 (which unlist handles by skipping)
        # To be safe, if nrow(Er) == 0, assign logical(0)
        if (nrow(Er) == 0) per_run_spike_TR_masks[[run_name]] <- logical(0)
        next
    }
    
    Er_t <- t(Er) # Voxels x Time_r
    
    # Lambda for this run's RPCA
    lambda_r <- if (current_opts$rpca_lambda_auto) {
      1 / sqrt(max(dim(Er_t))) 
    } else {
      if (is.null(current_opts$rpca_lambda_fixed)) stop("rpca_lambda_fixed must be provided if rpca_lambda_auto is FALSE.")
      current_opts$rpca_lambda_fixed
    }
    
    # Prepare arguments for rpca::rpca
    # We will only pass 'mu' if current_opts$rpca_mu is not NULL,
    # allowing rpca::rpca to use its internal default and auto-tuning otherwise.
    rpca_call_args <- list(
      M = Er_t,
      lambda = lambda_r,
      term.delta = current_opts$rpca_term_delta,
      max.iter = current_opts$rpca_max_iter,
      trace = current_opts$rpca_trace
    )
    
    if (!is.null(current_opts$rpca_mu)) {
      # Only add mu to the call if it's specified by the user
      # Otherwise, rpca package will use its own default: prod(dim(A))/(4*sum(abs(A))) and may auto-tune.
      # Explicitly calculating mu like: mu_r <- prod(dim(Er_t)) / (4 * sum(abs(Er_t)))
      # and passing it seemed to cause issues in tests, possibly overriding internal tuning.
      rpca_call_args$mu <- current_opts$rpca_mu
      message(sprintf("  Run %s: Using user-specified mu = %f for rpca.", run_name, current_opts$rpca_mu))
    } else {
      # Calculate our 'default' mu for logging/messaging if needed, but don't pass it to rpca()
      # This way, we rely on rpca's internal default mu handling.
      temp_mu_for_logging <- NA
      if (sum(abs(Er_t)) > 1e-9) {
        temp_mu_for_logging <- prod(dim(Er_t)) / (4 * sum(abs(Er_t)))
      } else {
        temp_mu_for_logging <- 1.0 # Fallback for logging if Er_t is zero
      }
      message(sprintf("  Run %s: rpca_mu is NULL, rpca::rpca will use its internal default mu (approx for logging: %.2e).", 
                      run_name, temp_mu_for_logging))
    }
    
    # k for this run's RPCA (target rank of L component of Er_t)
    # Should not exceed min(dim(Er_t))
    k_this_run <- min(current_opts$k_per_run_target, min(dim(Er_t)))
    if (k_this_run <= 0) {
        warning(sprintf("Cannot perform RPCA for run %s: k_this_run (%d) is not positive after adjustment. Skipping.", run_name, k_this_run))
        V_list[[run_name]] <- NULL
        glitch_ratios_per_run[run_name] <- NA
        next
    }
    
    message(sprintf("  Processing run %s (data: %d voxels x %d TRs, target k: %d, lambda: %.2e)", 
                    run_name, nrow(Er_t), ncol(Er_t), k_this_run, lambda_r))
    
    rpca_res_r <- NULL
    tryCatch({
      rpca_res_r <- do.call(rpca, rpca_call_args)
    }, error = function(e) {
      warning(sprintf("rpca::rpca failed for run %s: %s", run_name, e$message))
      rpca_res_r <<- NULL
    })

    if (is.null(rpca_res_r) || is.null(rpca_res_r$L)) {
      warning(sprintf("Robust PCA failed to produce L for run %s (or rpca_res_r is NULL).", run_name))
      V_list[[run_name]] <- NULL
      glitch_ratios_per_run[run_name] <- NA
      next
    }
    
    # As per proposal addendum: Vr <- rp$U (where rp = rpca(Er_t, ...))
    # This was based on a generic rpca call. The rpca package might return L and S,
    # and we need to SVD L to get the voxel-space components (U of L if L is Voxels x Time_r).
    # If L_r_t = rpca_res_r$L (Voxels x Time_r), then its left singular vectors are V_r.
    L_r_t <- rpca_res_r$L # Voxels x Time_r
    
    if (is.null(L_r_t) || min(dim(L_r_t)) == 0 || sum(abs(L_r_t)) < 1e-9) {
        warning(sprintf("RPCA for run %s yielded an empty or zero L component.", run_name))
        V_list[[run_name]] <- NULL
        glitch_ratios_per_run[run_name] <- NA
        next
    }

    # Perform SVD on L_r_t to get V_r (left singular vectors of L_r_t)
    # k_this_run is the target number of components from this L_r_t
    svd_L_r_t <- NULL
    k_eff_svd_L <- min(k_this_run, nrow(L_r_t), ncol(L_r_t))
    if (k_eff_svd_L < k_this_run) {
        message(sprintf("  Run %s: Effective k for SVD of L_r_t (%d) is less than k_this_run (%d) due to matrix dimensions.", 
                        run_name, k_eff_svd_L, k_this_run))
    }
    
    tryCatch({
        # nu = k_this_run means we want up to k_this_run left singular vectors
        svd_L_r_t <- svd(L_r_t, nu = k_eff_svd_L, nv = 0) 
    }, error = function(e) {
        warning(sprintf("SVD on L_r_t for run %s failed: %s", run_name, e$message))
        svd_L_r_t <<- NULL
    })
    
    if (is.null(svd_L_r_t) || is.null(svd_L_r_t$u) || ncol(svd_L_r_t$u) == 0) {
      warning(sprintf("SVD of L_r_t for run %s yielded no U components for V_r.", run_name))
      V_list[[run_name]] <- NULL
      glitch_ratios_per_run[run_name] <- NA
      next
    }
    V_r <- svd_L_r_t$u # Voxels x (up to) k_this_run
    
    # The original proposal note Vr <- rp$U might refer to an rpca implementation where $U 
    # directly gives the left singular vectors of L if L is the primary low-rank matrix of interest.
    # With rpca::rpca, we get L and S, then SVD L.
    
    if (is.null(V_r) || ncol(V_r) == 0) {
        warning(sprintf("RPCA for run %s yielded no voxel-space components (V_r is NULL or empty).", run_name))
        V_list[[run_name]] <- NULL
        glitch_ratios_per_run[run_name] <- NA
        next
    }
    
    V_list[[run_name]] <- V_r
    
    # Glitch ratio for this run
    S_r_t <- rpca_res_r$S
    L_r_t <- rpca_res_r$L
    energy_S_r <- sum(S_r_t^2)
    energy_L_r <- sum(L_r_t^2)
    if (energy_L_r > 1e-9) {
      glitch_ratios_per_run[run_name] <- energy_S_r / energy_L_r
    } else {
      glitch_ratios_per_run[run_name] <- NA
    }

    # Spike TR mask for this run using median |S| across voxels
    if (!is.null(S_r_t) && length(S_r_t) > 0 && sum(abs(S_r_t), na.rm = TRUE) > 1e-9) {
      # Only update if S_r_t is valid and non-zero
      if (sum(abs(S_r_t), na.rm = TRUE) < 1e-9) {
        # This condition is a bit redundant due to the outer if, but safe
        # spike_TR_mask_run <- rep(FALSE, ncol(S_r_t)) # ncol(S_r_t) is Time_r
        # per_run_spike_TR_masks[[run_name]] <- as.logical(spike_TR_mask_run) # Stays FALSE
      } else {
        s_t_star_run <- apply(abs(S_r_t), 2, stats::median, na.rm = TRUE)
        mad_s_t_star <- stats::mad(s_t_star_run, constant = 1, na.rm = TRUE)
        if (!is.finite(mad_s_t_star) || mad_s_t_star < 1e-9) {
          thresh_val <- stats::quantile(s_t_star_run,
                                        probs = current_opts$rpca_spike_percentile_thresh,
                                        na.rm = TRUE)
        } else {
          thresh_val <- median(s_t_star_run, na.rm = TRUE) +
            current_opts$rpca_spike_mad_thresh * mad_s_t_star
        }
        spike_TR_mask_run <- s_t_star_run > thresh_val
        if (length(spike_TR_mask_run) == nrow(Er)) { # Ensure correct length
             per_run_spike_TR_masks[[run_name]] <- as.logical(spike_TR_mask_run)
        } else {
            warning(sprintf("Run %s: Generated spike_TR_mask_run length (%d) did not match TRs for run (%d). Using all FALSEs.", 
                            run_name, length(spike_TR_mask_run), nrow(Er)))
            # per_run_spike_TR_masks[[run_name]] remains rep(FALSE, nrow(Er))
        }
      }
    } # else, it remains the pre-initialized rep(FALSE, nrow(Er))

    # After S_r_t <- rpca_res_r$S and L_r_t <- rpca_res_r$L
    if (!is.null(S_r_t)) {
        S_matrix_list_per_run_TpV[[run_name]] <- t(S_r_t) # Transpose S_r_t (Voxels x Time_r) to Time_r x Voxels
    } else {
        S_matrix_list_per_run_TpV[[run_name]] <- matrix(0, nrow = nrow(Er), ncol = ncol(Er)) # Placeholder if S is NULL
    }
  } # End per-run RPCA loop
  
  # Filter out NULLs from V_list (runs that failed RPCA)
  V_list_valid <- Filter(Negate(is.null), V_list)
  if (length(V_list_valid) == 0) {
    warning("RPCA failed for all runs, or no components were extracted. Cannot proceed.")
    return(NULL)
  }
  message(sprintf("Successfully extracted initial voxel components (V_r) from %d runs.", length(V_list_valid)))

  # --- 3. Merge Voxel-Space Components (V_r) to get V_global --- 
  V_global <- NULL
  k_actual_global_components <- 0

  if (current_opts$rpca_merge_strategy == "iterative") {
    message("Using Iterative Grassmann Averaging for V_global...")
    V_global <- .grassmann_merge_iterative(V_list_valid, k_global_target)
    if (!is.null(V_global)) {
      k_actual_global_components <- ncol(V_global)
    } else {
      warning("Iterative Grassmann Averaging failed to produce V_global.")
      return(NULL) # Or fallback to concat_svd if desired as a robust measure?
    }
  } else if (current_opts$rpca_merge_strategy == "concat_svd") {
    message("Using Concatenate & SVD method for V_global...")
    V_all_concat_list <- Filter(function(x) !is.null(x) && ncol(x) > 0, V_list_valid)
    if (length(V_all_concat_list) == 0) {
        warning("No valid V_r components to concatenate for SVD V_global step.")
        return(NULL)
    }
    V_all_concat <- do.call(cbind, V_all_concat_list)
    
    if (is.null(V_all_concat) || ncol(V_all_concat) == 0) {
        warning("Concatenation of V_r components for SVD V_global resulted in an empty matrix.")
        return(NULL)
    }
    
    message(sprintf("Performing SVD on concatenated V_r matrix (dims: %s) to find V_global (%d components)", 
                    paste(dim(V_all_concat), collapse="x"), k_global_target))
    
    k_for_svd_V_all <- min(k_global_target, ncol(V_all_concat), nrow(V_all_concat)) 
    if (k_for_svd_V_all <= 0) {
        warning(sprintf("Cannot perform SVD on concatenated V_r: k_for_svd_V_all (%d) is not positive.", k_for_svd_V_all))
        return(NULL)
    }
    
    svd_V_all <- NULL
    tryCatch({
        svd_V_all <- svd(V_all_concat, nu = k_for_svd_V_all, nv = 0)
    }, error = function(e){
        warning(paste("SVD on concatenated V_r components for V_global failed:", e$message))
        return(NULL)
    })
    
    if(is.null(svd_V_all) || is.null(svd_V_all$u) || ncol(svd_V_all$u) == 0) {
        warning("SVD on concatenated V_r for V_global failed to produce U vectors.")
        return(NULL)
    }
    V_global <- svd_V_all$u 
    k_actual_global_components <- ncol(V_global)
    message(sprintf("  V_global (concat_svd) obtained with %d components.", k_actual_global_components))

    # IF using concat_svd strategy, capture singular values:
    if (current_opts$rpca_merge_strategy == "concat_svd" && !is.null(svd_V_all)) {
      V_global_singular_values <- svd_V_all$d
    } else if (current_opts$rpca_merge_strategy == "iterative" && !is.null(V_global)) {
      # For iterative, a representative set of singular values would be harder to get simply.
      # One option: SVD the final V_global (which should be orthonormal) against itself to get values if meaningful,
      # or acknowledge that direct singular values for adaptation are best from concat_svd method.
      # For now, leave NULL if iterative, or compute SVD of final V_global if it makes sense.
      # Placeholder: if needed, svd(V_global)$d - but V_global cols are already principal components.
      # The singular values needed are from the matrix whose rank is being adapted.
      # This is usually the L component or the data matrix itself before RPCA.
      # The current Auto_Adapt_RPCA_Rank in workflow uses SVD of Y_residuals_current (placeholder).
      # For true adaptation, this function should output singular values from the L components before merging, or from merged V_all_concat.
      # For now, only populate from concat_svd strategy.
      if (verbose && current_opts$rpca_merge_strategy == "iterative") {
          message("Singular values for rank adaptation not directly available from iterative merge strategy in this function.")
      }
    }
  } else {
    stop(sprintf("Invalid rpca_merge_strategy: '%s'. Choose 'concat_svd' or 'iterative'.", current_opts$rpca_merge_strategy))
  }

  if (is.null(V_global) || k_actual_global_components == 0) {
      warning("V_global is NULL or has zero components after merging. Cannot proceed.")
      return(NULL)
  }

  # --- 4. Form Run-Specific Low-Rank Nuisance Time Courses (Cr) --- 
  C_r_list <- list()
  message("Forming run-specific temporal components C_r...")
  for (r_idx in seq_along(Y_residuals_list)) {
    run_name <- names(Y_residuals_list)[r_idx]
    Er <- Y_residuals_list[[r_idx]] 
    
    if (is.null(Er) || nrow(Er) == 0) {
        # Ensure C_r contributes zero rows but maintains correct ncol for rbind
        original_run_TRs <- ifelse(is.null(Y_residuals_list[[run_name]]), 0, nrow(Y_residuals_list[[run_name]]))
        C_r_list[[run_name]] <- matrix(0, nrow = original_run_TRs, ncol = k_actual_global_components) 
        next
    }
    
    C_r <- Er %*% V_global 
    C_r_list[[run_name]] <- C_r
  }
  
  C_components_cat <- do.call(rbind, C_r_list[paste0("run_", unique_runs)])
  
  if (is.null(C_components_cat) || nrow(C_components_cat) != nrow(Y_residuals_cat) || ncol(C_components_cat) != k_actual_global_components) {
      warning("Final concatenated C_components matrix dimensions are incorrect or matrix is NULL. Check C_r formation.")
      return(NULL)
  }
  
  # --- Concatenate S_matrix per run ---
  S_matrix_cat <- do.call(rbind, S_matrix_list_per_run_TpV[paste0("run_", unique_runs)])
  if (is.null(S_matrix_cat) || nrow(S_matrix_cat) != nrow(Y_residuals_cat) || ncol(S_matrix_cat) != ncol(Y_residuals_cat)) {
      warning("Concatenated S_matrix dimensions are incorrect or matrix is NULL. Returning NULL for S_matrix.")
      S_matrix_cat <- NULL # Ensure it's NULL if problematic
  }

  # Combine spike masks in original run order
  spike_TR_mask <- unlist(per_run_spike_TR_masks[paste0("run_", unique_runs)], use.names = FALSE)
  if (is.null(spike_TR_mask) || length(spike_TR_mask) != nrow(Y_residuals_cat)) {
    warning(sprintf("Global spike_TR_mask length (%d) mismatch with total timepoints (%d). Defaulting to all FALSE.", 
                    length(spike_TR_mask %||% 0), nrow(Y_residuals_cat)))
    spike_TR_mask <- rep(FALSE, nrow(Y_residuals_cat))
  } else {
    spike_TR_mask <- as.logical(spike_TR_mask)
  }

  message(sprintf("Multi-run RPCA: Returning %d concatenated temporal components.", ncol(C_components_cat %||% matrix(ncol=0))))
  return(list(C_components = C_components_cat,
              spike_TR_mask = spike_TR_mask,
              S_matrix_cat = S_matrix_cat,
              V_global_singular_values = V_global_singular_values
              ))
}

#' Iterative Grassmann Averaging of Voxel-Space Components
#'
#' Merges a list of voxel-space component matrices (V_r from each run/bag)
#' using an iterative Grassmann averaging approach to find a global V basis.
#'
#' @param V_list_valid A list of valid V_r matrices (Voxels x k_r). Each matrix should have the same number of rows (Voxels).
#' @param k_target_global Integer, the desired number of dimensions (columns) for the final V_global.
#' @return A matrix V_global (Voxels x k_target_global), or NULL if merging fails or k_target_global is 0.
#' @keywords internal
.grassmann_merge_iterative <- function(V_list_valid, k_target_global) {
  if (length(V_list_valid) == 0) {
    warning(".grassmann_merge_iterative: V_list_valid is empty.")
    return(NULL)
  }
  if (k_target_global <= 0) {
    warning(".grassmann_merge_iterative: k_target_global must be positive.")
    return(NULL)
  }

  # Initialize V_global with the first valid V_r, ensuring it has k_target_global or fewer columns
  V_global <- V_list_valid[[1]]
  if (ncol(V_global) > k_target_global) {
    V_global <- V_global[, 1:k_target_global, drop = FALSE]
  } else if (ncol(V_global) < k_target_global) {
    # This means the first run had fewer components than desired globally.
    # k_target_global might need to be adjusted, or this run might not be ideal for init.
    # For now, we proceed, the SVD later will be limited by ncol(M_proj)
    message(sprintf("Iterative merge init: First V_r has %d components, less than k_target_global %d.", ncol(V_global), k_target_global))
  }
  if (ncol(V_global) == 0) { # check if after truncation V_global is empty
      warning("Iterative merge init: V_global became empty after initial truncation/selection from first V_r.")
      return(NULL)
  }
  
  num_voxels <- nrow(V_global)

  if (length(V_list_valid) > 1) {
    for (r_idx in 2:length(V_list_valid)) {
      Vr <- V_list_valid[[r_idx]]
      if (is.null(Vr) || ncol(Vr) == 0 || nrow(Vr) != num_voxels) {
        warning(sprintf("Iterative merge: Skipping V_r for run/item %d due to NULL, 0 columns, or mismatched voxel dim (%d vs %d).", 
                        r_idx, ifelse(is.null(Vr), NA, nrow(Vr)), num_voxels))
        next
      }
      
      # Form a basis for the union of the two subspaces V_global and Vr
      # P_union will have at most rank(V_global) + rank(Vr) columns
      # Ensure V_global and Vr passed to cbind are not empty
      if (ncol(V_global) == 0) { # If V_global somehow became empty, re-initialize with current Vr
          V_global <- Vr
          if (ncol(V_global) > k_target_global) V_global <- V_global[, 1:k_target_global, drop=FALSE]
          if (ncol(V_global) == 0) { warning("Iterative merge: V_global re-init failed."); return(NULL); }
          next # Continue to next Vr if any
      }
      
      P_union_qr <- qr(cbind(V_global, Vr))
      # Check rank to avoid issues with qr.Q if matrix is all zeros or rank is tiny
      if (P_union_qr$rank == 0) {
          warning(sprintf("Iterative merge: cbind(V_global, Vr) for item %d is rank 0. Skipping merge step.", r_idx))
          next
      }
      P_union <- qr.Q(P_union_qr)
      
      if (ncol(P_union) == 0) {
           warning(sprintf("Iterative merge: P_union for item %d has 0 columns. Skipping merge step.", r_idx))
           next
      }

      # Average projection operator in the union basis
      M_proj <- t(P_union) %*% (V_global %*% t(V_global) + Vr %*% t(Vr)) %*% P_union
      
      # Number of components to extract from this merge step's SVD
      # Should be k_target_global, but not more than available from M_proj
      k_for_this_svd <- min(k_target_global, ncol(M_proj), nrow(M_proj)) 
      
      if (k_for_this_svd <= 0) {
        warning(sprintf("Iterative merge: k_for_this_svd is %d for item %d. Cannot perform SVD. Retaining previous V_global.", k_for_this_svd, r_idx))
        # V_global remains as it was before this problematic Vr
      } else {
        svd_M <- svd(M_proj, nu = k_for_this_svd, nv = 0)
        if (is.null(svd_M$u) || ncol(svd_M$u) == 0) {
            warning(sprintf("Iterative merge: SVD of M_proj for item %d yielded no U components. Retaining previous V_global.", r_idx))
        } else {
            V_global <- P_union %*% svd_M$u
        }
      }
    }
  }
  
  # Final check on V_global dimensions
  if (is.null(V_global) || ncol(V_global) == 0) {
      warning("Iterative merge: Final V_global is NULL or has 0 columns.")
      return(NULL)
  }
  
  # Ensure V_global has exactly k_target_global columns if possible, or fewer if rank was deficient
  if (ncol(V_global) > k_target_global) {
    V_global <- V_global[, 1:k_target_global, drop = FALSE]
  } else if (ncol(V_global) < k_target_global) {
    message(sprintf("Iterative merge: Final V_global has %d components, less than target %d.", ncol(V_global), k_target_global))
  }
  
  message(sprintf("Iterative Grassmann Merge: Final V_global obtained with %d components.", ncol(V_global)))
  return(V_global)
}

# Placeholder for the actual Grassmann iterative merge if needed later
# .grassmann_average_subspaces <- function(V_list, k_target) { ... }

#' Auto-Adaptive Rank Selection for RPCA-L
#'
#' Determines an appropriate rank for the low-rank component of RPCA
#' based on the "elbow" in a provided singular value spectrum.
#' The first index where the singular value ratio drops below
#' `drop_ratio` is chosen, then clamped between `k_min` and `k_max`.
#'
#' @param singular_values Numeric vector of singular values sorted in
#'   decreasing order.
#' @param drop_ratio Numeric scalar specifying the drop-off threshold
#'   relative to the first singular value. Default is `0.02` (2%).
#' @param k_min Integer minimum allowable rank. Default is `20`.
#' @param k_max Integer maximum allowable rank. Default is `50`.
#' @return Integer adaptive rank value within `[k_min, k_max]` and not
#'   exceeding `length(singular_values)`.
#' @export
Auto_Adapt_RPCA_Rank <- function(singular_values,
                                 drop_ratio = 0.02,
                                 k_min = 20L,
                                 k_max = 50L) {
  if (is.null(singular_values) || length(singular_values) == 0) {
    warning("Auto_Adapt_RPCA_Rank: singular_values is NULL or empty.")
    return(as.integer(k_min))
  }
  if (!is.numeric(singular_values)) {
    stop("singular_values must be numeric")
  }
  if (!is.finite(singular_values[1])) {
    warning("Auto_Adapt_RPCA_Rank: first singular value is not finite.")
    return(as.integer(k_min))
  }

  ratios <- singular_values / singular_values[1]
  k_candidate <- which(ratios < drop_ratio)[1]
  if (is.na(k_candidate)) {
    k_candidate <- length(singular_values)
  }

  k_candidate <- max(k_candidate, k_min)
  k_candidate <- min(k_candidate, k_max)
  k_candidate <- min(k_candidate, length(singular_values))

  return(as.integer(k_candidate))
}
