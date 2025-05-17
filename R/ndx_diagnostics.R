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

  # ---- Adaptive hyperparameters ----
  last_diag <- workflow_output$diagnostics_per_pass[[workflow_output$num_passes_completed]]

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

  html_lines <- c(html_lines, "<h2>Key Adaptive Hyperparameters</h2>", "<ul>")
  if (!is.null(last_diag$k_rpca_global)) {
    html_lines <- c(html_lines, sprintf("<li>k_rpca_global: %s</li>", last_diag$k_rpca_global))
  }
  if (!is.null(last_diag$num_spectral_sines)) {
    html_lines <- c(html_lines, sprintf("<li>num_spectral_sines: %s</li>", last_diag$num_spectral_sines))
  }
  if (!is.null(last_diag$lambda_parallel_noise)) {
    html_lines <- c(html_lines, sprintf("<li>lambda_parallel_noise: %.3f</li>", last_diag$lambda_parallel_noise))
  }
  if (!is.null(last_diag$lambda_perp_signal)) {
    html_lines <- c(html_lines, sprintf("<li>lambda_perp_signal: %.3f</li>", last_diag$lambda_perp_signal))
  }
  html_lines <- c(html_lines, "</ul>", "</body>", "</html>")

  html_path <- file.path(output_dir, filename)
  writeLines(html_lines, con = html_path)
  invisible(html_path)
}
