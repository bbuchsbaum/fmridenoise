#' Generate ND-X HTML Diagnostic Report
#'
#' This function creates a simple HTML report summarizing key diagnostics from
#' a run of `NDX_Process_Subject`. It plots the Denoising Efficacy Score (DES)
#' across passes, compares the power spectral density (PSD) of residuals from
#' Pass 0 and the final ND-X pass, and optionally shows a spike "carpet" plot
#' from the RPCA `S` matrix.
#'
#' @param workflow_output List returned by `NDX_Process_Subject`.
#' @param pass0_residuals Matrix of residuals from `ndx_initial_glm` used as the
#'   starting point for ND-X. Required for the PSD comparison.
#' @param TR Numeric repetition time (seconds) for frequency axis in the PSD plot.
#' @param output_dir Directory in which to save the HTML file and PNG figures.
#' @param filename Name of the HTML file to create. Default "ndx_diagnostic_report.html".
#' @return Invisibly, the path to the generated HTML file.
#' @export ndx_generate_html_report
ndx_generate_html_report <- function(workflow_output,
                                     pass0_residuals,
                                     TR,
                                     output_dir = "./diagnostics",
                                     filename = "ndx_diagnostic_report.html") {
  if (is.null(workflow_output) || !is.list(workflow_output)) {
    stop("workflow_output must be a list from NDX_Process_Subject")
  }
  if (!is.matrix(pass0_residuals)) {
    stop("pass0_residuals must be a numeric matrix")
  }
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a positive numeric scalar")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---- DES per pass plot ----
  des_values <- sapply(workflow_output$diagnostics_per_pass, function(d) d$DES)
  des_png <- file.path(output_dir, "des_per_pass.png")
  grDevices::png(des_png, width = 600, height = 400)
  graphics::plot(seq_along(des_values), des_values, type = "b",
                 xlab = "Pass", ylab = "DES",
                 main = "Denoising Efficacy Score per Pass")
  grDevices::dev.off()

  # ---- Residual PSD plot ----
  psd_png <- file.path(output_dir, "residual_psd.png")
  if (!is.null(workflow_output$Y_residuals_final_unwhitened)) {
    res0_mean <- rowMeans(pass0_residuals, na.rm = TRUE)
    res_final_mean <- rowMeans(workflow_output$Y_residuals_final_unwhitened, na.rm = TRUE)
    psd0 <- psd::pspectrum(res0_mean, x.frqsamp = 1/TR, plot = FALSE)
    psdF <- psd::pspectrum(res_final_mean, x.frqsamp = 1/TR, plot = FALSE)
    grDevices::png(psd_png, width = 600, height = 400)
    graphics::plot(psd0$freq, psd0$spec, type = "l", col = "red",
                   xlab = "Frequency (Hz)", ylab = "PSD",
                   main = "Residual Power Spectral Density", log = "y")
    graphics::lines(psdF$freq, psdF$spec, col = "blue")
    graphics::legend("topright", legend = c("Pass 0", "Final"),
                     col = c("red", "blue"), lty = 1, bty = "n")
    grDevices::dev.off()
  } else {
    psd_png <- NULL
  }

  # ---- Spike carpet plot ----
  carpet_png <- NULL
  if (!is.null(workflow_output$S_matrix_rpca_final)) {
    carpet_png <- file.path(output_dir, "spike_carpet.png")
    S_mat <- workflow_output$S_matrix_rpca_final
    grDevices::png(carpet_png, width = 600, height = 400)
    graphics::image(t(abs(S_mat) > 0), axes = FALSE, useRaster = TRUE,
                    xlab = "Time", ylab = "Voxel",
                    main = "RPCA Spike Carpet")
    grDevices::box()
    grDevices::dev.off()
  }

  # ---- Annihilation Mode diagnostics plot (if enabled) ----
  gdlite_png <- NULL
  if (!is.null(workflow_output$annihilation_mode_active) && workflow_output$annihilation_mode_active &&
      !is.null(workflow_output$gdlite_diagnostics) && !is.null(workflow_output$gdlite_diagnostics$r2_vals_by_k)) {
    
    gdlite_png <- file.path(output_dir, "gdlite_r2_by_k.png")
    r2_vals <- workflow_output$gdlite_diagnostics$r2_vals_by_k
    optimal_k <- workflow_output$gdlite_diagnostics$optimal_k
    k_vals <- seq_along(r2_vals)
    
    grDevices::png(gdlite_png, width = 600, height = 400)
    graphics::plot(k_vals, r2_vals, type = "b", col = "darkgreen",
                  xlab = "Number of PCs (k)", ylab = "Cross-validated R²",
                  main = "GLMdenoise-Lite: Cross-validated R² by Number of PCs")
    
    # Highlight optimal k if available
    if (!is.null(optimal_k) && optimal_k > 0 && optimal_k <= length(r2_vals)) {
      graphics::points(optimal_k, r2_vals[optimal_k], pch = 19, col = "red", cex = 1.5)
      graphics::text(optimal_k, r2_vals[optimal_k], 
                    labels = sprintf("k=%d", optimal_k),
                    pos = 3, col = "red")
    }
    
    grDevices::dev.off()
  }

  # ---- Adaptive hyperparameters ----
  # Get the diagnostics from the last completed pass
  last_pass_num <- workflow_output$num_passes_completed
  last_diag <- workflow_output$diagnostics_per_pass[[last_pass_num]]

  html_lines <- c(
    "<html>",
    "<head><title>ND-X Diagnostic Report</title></head>",
    "<body>",
    "<h1>ND-X Diagnostic Report</h1>",
    "<h2>Denoising Efficacy Score per Pass</h2>",
    sprintf("<img src='%s' alt='DES per pass'>", basename(des_png))
  )

  if (!is.null(psd_png)) {
    html_lines <- c(html_lines,
      "<h2>Residual Power Spectral Density</h2>",
      sprintf("<img src='%s' alt='Residual PSD'>", basename(psd_png)))
  }

  if (!is.null(carpet_png)) {
    html_lines <- c(html_lines,
      "<h2>Spike Carpet Plot</h2>",
      sprintf("<img src='%s' alt='Spike carpet'>", basename(carpet_png)))
  }

  # Add Annihilation Mode diagnostics if active
  if (!is.null(workflow_output$annihilation_mode_active) && workflow_output$annihilation_mode_active) {
    html_lines <- c(html_lines, "<h2>Annihilation Mode Diagnostics</h2>")
    
    if (!is.null(gdlite_png)) {
      html_lines <- c(html_lines, 
        sprintf("<img src='%s' alt='GLMdenoise-Lite R² by k'>", basename(gdlite_png)))
    }
    
    html_lines <- c(html_lines, "<ul>")
    
    if (!is.null(workflow_output$gdlite_pcs)) {
      html_lines <- c(html_lines, 
        sprintf("<li>Number of GLMdenoise PCs used: %d</li>", ncol(workflow_output$gdlite_pcs)))
    }
    
    if (!is.null(workflow_output$gdlite_diagnostics)) {
      if (!is.null(workflow_output$gdlite_diagnostics$optimal_k)) {
        html_lines <- c(html_lines, 
          sprintf("<li>Optimal number of GLMdenoise PCs: %d</li>", workflow_output$gdlite_diagnostics$optimal_k))
      }
      
      if (!is.null(workflow_output$gdlite_diagnostics$noise_pool_mask)) {
        noise_pool_size <- sum(workflow_output$gdlite_diagnostics$noise_pool_mask)
        total_voxels <- length(workflow_output$gdlite_diagnostics$noise_pool_mask)
        html_lines <- c(html_lines, 
          sprintf("<li>Noise pool voxels: %d (%.1f%% of total)</li>", 
                  noise_pool_size, 100 * noise_pool_size / total_voxels))
      }
      
      if (!is.null(workflow_output$gdlite_diagnostics$good_voxels_mask)) {
        good_voxels_size <- sum(workflow_output$gdlite_diagnostics$good_voxels_mask)
        total_voxels <- length(workflow_output$gdlite_diagnostics$good_voxels_mask)
        html_lines <- c(html_lines, 
          sprintf("<li>Good voxels: %d (%.1f%% of total)</li>", 
                  good_voxels_size, 100 * good_voxels_size / total_voxels))
      }
    }
    
    html_lines <- c(html_lines, "</ul>")
  }

  html_lines <- c(html_lines, "<h2>Key Adaptive Hyperparameters</h2>", "<ul>")
  
  # Extract hyperparameters from the correct location in last_diag
  if (!is.null(last_diag$pass_options) && !is.null(last_diag$pass_options$current_rpca_k_target)) {
    html_lines <- c(html_lines, 
      sprintf("<li>k_rpca_global: %d</li>", last_diag$pass_options$current_rpca_k_target))
  } else if (!is.null(last_diag$k_rpca_global)) {
    html_lines <- c(html_lines, 
      sprintf("<li>k_rpca_global: %d</li>", last_diag$k_rpca_global))
  }
  
  if (!is.null(last_diag$num_spectral_sines)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Number of spectral sine pairs: %d</li>", last_diag$num_spectral_sines))
  }
  
  # Add ridge regression lambda parameters
  if (!is.null(last_diag$lambda_parallel_noise)) {
    html_lines <- c(html_lines, 
      sprintf("<li>lambda_parallel_noise (noise regressors): %.3f</li>", last_diag$lambda_parallel_noise))
  }
  
  if (!is.null(last_diag$lambda_perp_signal)) {
    html_lines <- c(html_lines, 
      sprintf("<li>lambda_perp_signal (task/signal regressors): %.3f</li>", last_diag$lambda_perp_signal))
  }
  
  # For Annihilation Mode, add unique component lambda if available
  if (!is.null(workflow_output$annihilation_mode_active) && 
      workflow_output$annihilation_mode_active &&
      !is.null(last_diag$pass_options$annihilation_active_this_pass) &&
      last_diag$pass_options$annihilation_active_this_pass) {
      
    html_lines <- c(html_lines, "<li><b>Annihilation Mode was active for this pass</b></li>")
  }
  
  # Add convergence info
  html_lines <- c(html_lines, 
    sprintf("<li>Number of passes completed: %d</li>", workflow_output$num_passes_completed))
  
  if (workflow_output$num_passes_completed < length(workflow_output$diagnostics_per_pass)) {
    html_lines <- c(html_lines, "<li>Convergence: Reached before maximum passes</li>")
  } else {
    html_lines <- c(html_lines, "<li>Convergence: Reached maximum passes</li>")
  }
  
  # Add final DES value
  if (!is.null(last_diag$DES)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Final Denoising Efficacy Score (DES): %.4f</li>", last_diag$DES))
  }
  
  # Add final rho noise projection if available
  if (!is.null(last_diag$rho_noise_projection)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Final Rho Noise Projection: %.4f</li>", last_diag$rho_noise_projection))
  }
  
  html_lines <- c(html_lines, "</ul>", "</body>", "</html>")

  html_path <- file.path(output_dir, filename)
  writeLines(html_lines, con = html_path)
  invisible(html_path)
}

#' Generate "ND-X Certified Clean" JSON Sidecar
#'
#' This function creates a JSON file summarizing key diagnostic metrics and
#' adaptive hyperparameters from an ND-X run. The JSON structure follows the
#' specification described in the ND-X proposal and is intended to be consumed
#' by downstream tools.
#'
#' @param workflow_output List returned by `NDX_Process_Subject`.
#' @param output_path File path to save the JSON sidecar. Defaults to
#'   "sub-ndx.json" in the current working directory.
#' @return Invisibly, the path to the written JSON file.
#' @export ndx_generate_json_certificate
ndx_generate_json_certificate <- function(workflow_output,
                                          output_path = "sub-ndx.json") {
  if (is.null(workflow_output) || !is.list(workflow_output)) {
    stop("workflow_output must be a list from NDX_Process_Subject")
  }
  if (is.null(workflow_output$diagnostics_per_pass) ||
      length(workflow_output$diagnostics_per_pass) == 0) {
    stop("workflow_output$diagnostics_per_pass is missing")
  }

  last_diag <- workflow_output$diagnostics_per_pass[[workflow_output$num_passes_completed]]

  val_or_na <- function(x) if (is.null(x)) NA else x

  cert_list <- list(
    ndx_version = as.character(utils::packageVersion("ndx")),
    final_DES = val_or_na(last_diag$DES),
    num_passes_converged = val_or_na(workflow_output$num_passes_completed),
    final_rho_noise_projection = val_or_na(last_diag$rho_noise_projection),
    final_adaptive_hyperparameters = list(
      k_rpca_global = val_or_na(last_diag$k_rpca_global),
      num_spectral_sines = val_or_na(last_diag$num_spectral_sines),
      lambda_parallel_noise = val_or_na(last_diag$lambda_parallel_noise),
      lambda_perp_signal = val_or_na(last_diag$lambda_perp_signal)
    )
  )

  jsonlite::write_json(cert_list, output_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(output_path)
}
