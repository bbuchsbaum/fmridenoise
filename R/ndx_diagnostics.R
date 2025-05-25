#' Generate ND-X HTML Diagnostic Report
#'
#' This function creates a simple HTML report summarizing key diagnostics from
#' a run of `NDX_Process_Subject`. It plots the Denoising Efficacy Score (DES)
#' across passes, compares the power spectral density (PSD) of residuals from
#' Pass 0 and the final ND-X pass, visualizes \eqn{\beta}-stability across runs,
#' shows the Ljung--Box p-value progression, and optionally displays a spike
#' "carpet" plot from the RPCA `S` matrix.
#'
#' @param workflow_output List returned by `NDX_Process_Subject`.
#' @param pass0_residuals Matrix of residuals from `ndx_initial_glm` used as the
#'   starting point for ND-X. Required for the PSD comparison.
#' @param TR Numeric repetition time (seconds) for frequency axis in the PSD plot.
#' @param output_dir Directory in which to save the HTML file and PNG figures.
#' @param filename Name of the HTML file to create. Default "ndx_diagnostic_report.html".
#' @param include_progressive_enhancement Logical, whether to include the four-panel
#'   progressive enhancement visualization (requires Annihilation Mode). Default TRUE.
#' @return Invisibly, the path to the generated HTML file.
#' @export ndx_generate_html_report
ndx_generate_html_report <- function(workflow_output,
                                     pass0_residuals,
                                     TR,
                                     output_dir = "./diagnostics",
                                     filename = "ndx_diagnostic_report.html",
                                     include_progressive_enhancement = TRUE) {
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

  # ---- Beta stability plot ----
  beta_png <- NULL
  if (!is.null(workflow_output$beta_history_per_pass)) {
    beta_vals <- tryCatch(
      calculate_beta_stability(workflow_output$beta_history_per_pass),
      error = function(e) rep(NA_real_, length(workflow_output$beta_history_per_pass))
    )
    beta_png <- file.path(output_dir, "beta_stability_per_pass.png")
    grDevices::png(beta_png, width = 600, height = 400)
    graphics::plot(seq_along(beta_vals), beta_vals, type = "b",
                   xlab = "Pass", ylab = "Beta Stability",
                   ylim = c(-1, 1),
                   main = "\u03B2-Stability Across Runs")
    grDevices::dev.off()
  }

  # ---- Residual PSD plot ----
  psd_png <- file.path(output_dir, "residual_psd.png")
  if (!is.null(workflow_output$Y_residuals_final_unwhitened)) {
    if (nrow(pass0_residuals) != nrow(workflow_output$Y_residuals_final_unwhitened)) {
      stop(
        sprintf(
          "Mismatch in residual lengths: pass0_residuals has %d rows, \
workflow_output$Y_residuals_final_unwhitened has %d rows",
          nrow(pass0_residuals),
          nrow(workflow_output$Y_residuals_final_unwhitened)
        )
      )
    }
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

  # ---- Ljung-Box p-value plot ----
  ljung_png <- NULL
  if (!is.null(workflow_output$diagnostics_per_pass)) {
    ljung_vals <- sapply(workflow_output$diagnostics_per_pass,
                         function(d) d$ljung_box_p)
    if (any(!is.na(ljung_vals))) {
      ljung_png <- file.path(output_dir, "ljung_box_pvalues.png")
      grDevices::png(ljung_png, width = 600, height = 400)
      graphics::plot(seq_along(ljung_vals), ljung_vals, type = "b",
                     xlab = "Pass", ylab = "Ljung-Box p-value",
                     ylim = c(0, 1),
                     main = "Residual Whiteness (Ljung-Box)")
      graphics::abline(h = 0.05, col = "red", lty = 2)
      grDevices::dev.off()
    }
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
    graphics::box()
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
    "<head>",
    "<title>ND-X Diagnostic Report</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }",
    "h1, h2, h3 { color: #333; }",
    ".section { margin: 20px 0; }",
    ".card { border: 1px solid #ddd; border-radius: 5px; padding: 15px; margin-bottom: 20px; }",
    ".params-list { list-style-type: none; padding: 0; }",
    ".params-list li { padding: 5px 0; border-bottom: 1px solid #eee; }",
    ".verdict { display: inline-block; padding: 10px 15px; border-radius: 5px; font-weight: bold; color: white; }",
    "</style>",
    "</head>",
    "<body>",
    "<h1>ND-X Diagnostic Report</h1>"
  )

  # ---- Add DES section ----
  html_lines <- c(html_lines,
    "<div class='section'>",
    "<h2>Denoising Efficacy Score per Pass</h2>",
    sprintf("<img src='%s' alt='DES per pass'>", basename(des_png)),
    "</div>"
  )

  # ---- Add beta stability section ----
  if (!is.null(beta_png)) {
    html_lines <- c(html_lines,
      "<div class='section'>",
      "<h2>&beta;-Stability Across Runs</h2>",
      sprintf("<img src='%s' alt='Beta stability'>", basename(beta_png)),
      "</div>"
    )
  }

  # ---- Add PSD section ----
  if (!is.null(psd_png)) {
    html_lines <- c(html_lines,
      "<div class='section'>",
      "<h2>Residual Power Spectral Density</h2>",
      sprintf("<img src='%s' alt='Residual PSD'>", basename(psd_png)),
      "</div>"
    )
  }

  # ---- Add Ljung-Box section ----
  if (!is.null(ljung_png)) {
    html_lines <- c(html_lines,
      "<div class='section'>",
      "<h2>Ljung-Box p-value per Pass</h2>",
      sprintf("<img src='%s' alt='Ljung-Box p-values'>", basename(ljung_png)),
      "</div>"
    )
  }

  # ---- Add Spike carpet section ----
  if (!is.null(carpet_png)) {
    html_lines <- c(html_lines,
      "<div class='section'>",
      "<h2>Spike Carpet Plot</h2>",
      sprintf("<img src='%s' alt='Spike carpet'>", basename(carpet_png)),
      "</div>"
    )
  }

  # ---- Add Annihilation Mode diagnostics if active ----
  if (!is.null(workflow_output$annihilation_mode_active) && workflow_output$annihilation_mode_active) {
    html_lines <- c(html_lines, "<div class='section'>", "<h2>Annihilation Mode Diagnostics</h2>")
    
    # Add GLMdenoise-Lite R² by k plot if available
    if (!is.null(gdlite_png)) {
      html_lines <- c(html_lines, 
        sprintf("<img src='%s' alt='GLMdenoise-Lite R² by k'>", basename(gdlite_png))
      )
    }
    
    # Add Annihilation Verdict
    if (!is.null(workflow_output$gdlite_pcs) && 
        (!is.null(workflow_output$rpca_orthogonalized) || !is.null(workflow_output$spectral_orthogonalized))) {
      
      # Calculate component variances
      var_gdlite <- sum(workflow_output$gdlite_pcs^2, na.rm = TRUE)
      
      var_ndx_unique <- 0
      if (!is.null(workflow_output$rpca_orthogonalized)) {
        var_ndx_unique <- var_ndx_unique + sum(workflow_output$rpca_orthogonalized^2, na.rm = TRUE)
      }
      if (!is.null(workflow_output$spectral_orthogonalized)) {
        var_ndx_unique <- var_ndx_unique + sum(workflow_output$spectral_orthogonalized^2, na.rm = TRUE)
      }
      
      verdict_ratio <- if (var_gdlite > 0) var_ndx_unique / var_gdlite else 0
      
      # Determine verdict category
      if (verdict_ratio < 0.1) {
        verdict <- "Tie"
        color <- "#888888"
        description <- "GLMdenoise and ND-X unique components capture similar noise variance."
      } else if (verdict_ratio < 0.5) {
        verdict <- "Win"
        color <- "#5cb85c"
        description <- "ND-X unique components capture additional noise not found by GLMdenoise."
      } else if (verdict_ratio < 1.0) {
        verdict <- "Decisive Win"
        color <- "#5bc0de"
        description <- "ND-X unique components capture substantial additional noise."
      } else {
        verdict <- "Annihilation"
        color <- "#d9534f"
        description <- "ND-X unique components capture more noise variance than GLMdenoise components."
      }
      
      html_lines <- c(html_lines,
        "<div class='card' style='text-align: center;'>",
        "<h3>Annihilation Verdict</h3>",
        sprintf("<div class='verdict' style='background-color: %s;'>%s</div>", color, verdict),
        sprintf("<p>Variance Ratio (ND-X Unique / GLMdenoise): <strong>%.2f</strong></p>", verdict_ratio),
        sprintf("<p style='font-style: italic;'>%s</p>", description),
        "</div>"
      )
    }
    
    # Add GLMdenoise stats
    html_lines <- c(html_lines, "<div class='card'>", "<h3>GLMdenoise-Lite Statistics</h3>", "<ul class='params-list'>")
    
    if (!is.null(workflow_output$gdlite_pcs)) {
      html_lines <- c(html_lines, 
        sprintf("<li>Number of GLMdenoise PCs used: <strong>%d</strong></li>", ncol(workflow_output$gdlite_pcs))
      )
    }
    
    if (!is.null(workflow_output$gdlite_diagnostics)) {
      if (!is.null(workflow_output$gdlite_diagnostics$optimal_k)) {
        html_lines <- c(html_lines, 
          sprintf("<li>Optimal number of GLMdenoise PCs: <strong>%d</strong></li>", 
                  workflow_output$gdlite_diagnostics$optimal_k)
        )
      }
      
      if (!is.null(workflow_output$gdlite_diagnostics$noise_pool_mask)) {
        noise_pool_size <- sum(workflow_output$gdlite_diagnostics$noise_pool_mask)
        total_voxels <- length(workflow_output$gdlite_diagnostics$noise_pool_mask)
        html_lines <- c(html_lines, 
          sprintf("<li>Noise pool voxels: <strong>%d</strong> (%.1f%% of total)</li>", 
                  noise_pool_size, 100 * noise_pool_size / total_voxels)
        )
      }
      
      if (!is.null(workflow_output$gdlite_diagnostics$good_voxels_mask)) {
        good_voxels_size <- sum(workflow_output$gdlite_diagnostics$good_voxels_mask)
        total_voxels <- length(workflow_output$gdlite_diagnostics$good_voxels_mask)
        html_lines <- c(html_lines, 
          sprintf("<li>Good voxels: <strong>%d</strong> (%.1f%% of total)</li>", 
                  good_voxels_size, 100 * good_voxels_size / total_voxels)
        )
      }
    }
    
    html_lines <- c(html_lines, "</ul>", "</div>", "</div>")
  }

  # ---- Add key adaptive hyperparameters section ----
  html_lines <- c(html_lines, 
    "<div class='section'>",
    "<h2>Key Adaptive Hyperparameters</h2>",
    "<div class='card'>",
    "<ul class='params-list'>"
  )
  
  # Extract hyperparameters from the correct location in last_diag
  if (!is.null(last_diag$pass_options) && !is.null(last_diag$pass_options$current_rpca_k_target)) {
    html_lines <- c(html_lines, 
      sprintf("<li>k_rpca_global: <strong>%d</strong></li>", last_diag$pass_options$current_rpca_k_target)
    )
  } else if (!is.null(last_diag$k_rpca_global)) {
    html_lines <- c(html_lines, 
      sprintf("<li>k_rpca_global: <strong>%d</strong></li>", last_diag$k_rpca_global)
    )
  }
  
  if (!is.null(last_diag$num_spectral_sines)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Number of spectral sine pairs: <strong>%d</strong></li>", last_diag$num_spectral_sines)
    )
  }
  
  # Add ridge regression lambda parameters
  if (!is.null(last_diag$lambda_parallel_noise)) {
    html_lines <- c(html_lines, 
      sprintf("<li>lambda_parallel_noise (noise regressors): <strong>%.3f</strong></li>", last_diag$lambda_parallel_noise)
    )
  }
  
  if (!is.null(last_diag$lambda_perp_signal)) {
    html_lines <- c(html_lines, 
      sprintf("<li>lambda_perp_signal (task/signal regressors): <strong>%.3f</strong></li>", last_diag$lambda_perp_signal)
    )
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
    sprintf("<li>Number of passes completed: <strong>%d</strong></li>", workflow_output$num_passes_completed)
  )
  
  if (workflow_output$num_passes_completed < length(workflow_output$diagnostics_per_pass)) {
    html_lines <- c(html_lines, "<li>Convergence: <strong>Reached before maximum passes</strong></li>")
  } else {
    html_lines <- c(html_lines, "<li>Convergence: <strong>Reached maximum passes</strong></li>")
  }
  
  # Add final DES value
  if (!is.null(last_diag$DES)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Final Denoising Efficacy Score (DES): <strong>%.4f</strong></li>", last_diag$DES)
    )
  }
  
  # Add final rho noise projection if available
  if (!is.null(last_diag$rho_noise_projection)) {
    html_lines <- c(html_lines, 
      sprintf("<li>Final Rho Noise Projection: <strong>%.4f</strong></li>", last_diag$rho_noise_projection)
    )
  }
  
  html_lines <- c(html_lines, "</ul>", "</div>", "</div>")

  # ---- Add Progressive Enhancement section if requested ----
  if (include_progressive_enhancement && 
      !is.null(workflow_output$annihilation_mode_active) && 
      workflow_output$annihilation_mode_active) {
    
    # Use the ndx_progressive_viz functions to generate this section
    progressive_html <- ndx_generate_progressive_enhancement(
      workflow_output = workflow_output,
      output_dir = output_dir
    )
    
    html_lines <- c(html_lines,
      "<div class='section'>",
      "<h2>Progressive Enhancement Visualization</h2>",
      "<p>This interactive visualization shows how each stage of ND-X processing improves denoising effectiveness.</p>",
      progressive_html,
      "</div>"
    )
  }

  # Close HTML document
  html_lines <- c(html_lines, "</body>", "</html>")

  # Write HTML file
  html_path <- file.path(output_dir, filename)
  writeLines(html_lines, con = html_path)
  invisible(html_path)
}

#' Generate "ND-X Certified Clean" JSON Sidecar
#'
#' This function creates a JSON file summarizing key diagnostic metrics and
#' adaptive hyperparameters from an ND-X run. The resulting JSON includes
#' fields for the ND-X version, final Denoising Efficacy Score (DES),
#' final beta stability, Ljung-Box p-value (from final pass), median Ljung-Box
#' p-value (from final unwhitened residuals), the Annihilation verdict
#' ratio and category (if applicable), convergence information (number of passes),
#' final rho noise projection, and all final adaptive hyperparameters.
#' Downstream tools can parse this certificate to verify denoising quality.
#' The workflow hash included in the JSON is computed from the diagnostics and
#' key settings rather than the full `workflow_output` object to avoid
#' excessive memory usage.
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

  verdict_stats <- ndx_annihilation_verdict_stats(workflow_output)

  cert_list <- list(
    ndx_version = as.character(utils::packageVersion("ndx")),
    final_DES = val_or_na(last_diag$DES),
    ljung_box_p = val_or_na(last_diag$ljung_box_p),
    var_ratio = val_or_na(verdict_stats$var_ratio),
    verdict = val_or_na(verdict_stats$verdict),
    passes_converged = val_or_na(workflow_output$num_passes_completed),
    final_rho_noise_projection = val_or_na(last_diag$rho_noise_projection),
    final_beta_stability = val_or_na(tail(calculate_beta_stability(workflow_output$beta_history_per_pass),1)),
    ljung_box_p_median = val_or_na(median(compute_ljung_box_pvalues(workflow_output$Y_residuals_final_unwhitened), na.rm=TRUE)),
    final_adaptive_hyperparameters = list(
      k_rpca_global = val_or_na(last_diag$k_rpca_global),
      num_hrf_clusters = val_or_na(last_diag$num_hrf_clusters),
      num_spectral_sines = val_or_na(last_diag$num_spectral_sines),
      lambda_parallel_noise = val_or_na(last_diag$lambda_parallel_noise),
      lambda_perp_signal = val_or_na(last_diag$lambda_perp_signal)
    ),
    timestamp = as.character(Sys.time()),
    workflow_hash = digest::digest(list(
      diagnostics = workflow_output$diagnostics_per_pass,
      settings = list(
        num_passes_completed = workflow_output$num_passes_completed,
        annihilation_mode_active = val_or_na(workflow_output$annihilation_mode_active)
      )
    ))
  )

  jsonlite::write_json(cert_list, output_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(output_path)
}
