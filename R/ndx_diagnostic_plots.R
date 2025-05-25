# Helper functions for generating diagnostic plots as PNG files
# These functions are internal and return the path to the generated file
# or NULL if the plot could not be produced.

#' Plot DES per pass and save to PNG
#' @keywords internal
ndx_plot_des_per_pass <- function(diagnostics_per_pass, output_dir) {
  if (!is.list(diagnostics_per_pass) || length(diagnostics_per_pass) == 0) {
    stop("diagnostics_per_pass must be a non-empty list")
  }
  des_vals <- vapply(diagnostics_per_pass, function(d) d$DES %||% NA_real_, numeric(1))
  png_file <- file.path(output_dir, "des_per_pass.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  graphics::plot(seq_along(des_vals), des_vals, type = "b", xlab = "Pass", ylab = "DES",
                 main = "Denoising Efficacy Score per Pass")
  grDevices::dev.off()
  png_file
}

#' Plot beta stability across passes
#' @keywords internal
ndx_plot_beta_stability <- function(beta_history_per_pass, output_dir) {
  if (!is.list(beta_history_per_pass) || length(beta_history_per_pass) == 0) {
    stop("beta_history_per_pass must be a non-empty list")
  }
  beta_vals <- calculate_beta_stability(beta_history_per_pass)
  png_file <- file.path(output_dir, "beta_stability_per_pass.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  graphics::plot(seq_along(beta_vals), beta_vals, type = "b", xlab = "Pass",
                 ylab = "Beta Stability", ylim = c(-1, 1),
                 main = "\u03B2-Stability Across Runs")
  grDevices::dev.off()
  png_file
}

#' Plot residual PSD comparison
#' @keywords internal
ndx_plot_residual_psd <- function(pass0_residuals, final_residuals, TR, output_dir) {
  if (!is.matrix(pass0_residuals) || !is.matrix(final_residuals)) {
    stop("pass0_residuals and final_residuals must be matrices")
  }
  if (nrow(pass0_residuals) != nrow(final_residuals)) {
    stop("Residual matrices must have the same number of rows")
  }
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a positive numeric scalar")
  }
  res0_mean <- rowMeans(pass0_residuals, na.rm = TRUE)
  res_final_mean <- rowMeans(final_residuals, na.rm = TRUE)
  psd0 <- psd::pspectrum(res0_mean, x.frqsamp = 1/TR, plot = FALSE)
  psdF <- psd::pspectrum(res_final_mean, x.frqsamp = 1/TR, plot = FALSE)
  png_file <- file.path(output_dir, "residual_psd.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  graphics::plot(psd0$freq, psd0$spec, type = "l", col = "red", log = "y",
                 xlab = "Frequency (Hz)", ylab = "PSD",
                 main = "Residual Power Spectral Density")
  graphics::lines(psdF$freq, psdF$spec, col = "blue")
  graphics::legend("topright", legend = c("Pass 0", "Final"), col = c("red", "blue"), lty = 1, bty = "n")
  grDevices::dev.off()
  png_file
}

#' Plot Ljung-Box p-values per pass
#' @keywords internal
ndx_plot_ljung_box_pvalues <- function(diagnostics_per_pass, output_dir) {
  if (!is.list(diagnostics_per_pass) || length(diagnostics_per_pass) == 0) {
    stop("diagnostics_per_pass must be a non-empty list")
  }
  ljung_vals <- vapply(diagnostics_per_pass, function(d) d$ljung_box_p %||% NA_real_, numeric(1))
  if (all(is.na(ljung_vals))) return(NULL)
  png_file <- file.path(output_dir, "ljung_box_pvalues.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  graphics::plot(seq_along(ljung_vals), ljung_vals, type = "b", xlab = "Pass",
                 ylab = "Ljung-Box p-value", ylim = c(0, 1),
                 main = "Residual Whiteness (Ljung-Box)")
  graphics::abline(h = 0.05, col = "red", lty = 2)
  grDevices::dev.off()
  png_file
}

#' Plot spike carpet from RPCA S matrix
#' @keywords internal
ndx_plot_spike_carpet <- function(S_matrix, output_dir) {
  if (!is.matrix(S_matrix)) {
    stop("S_matrix must be a matrix")
  }
  png_file <- file.path(output_dir, "spike_carpet.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  graphics::image(t(abs(S_matrix) > 0), axes = FALSE, useRaster = TRUE,
                  xlab = "Time", ylab = "Voxel", main = "RPCA Spike Carpet")
  graphics::box()
  grDevices::dev.off()
  png_file
}

#' Plot GLMdenoise R2 by k
#' @keywords internal
ndx_plot_gdlite_r2_by_k <- function(gdlite_diagnostics, output_dir) {
  if (is.null(gdlite_diagnostics$r2_vals_by_k)) {
    stop("gdlite_diagnostics must contain r2_vals_by_k")
  }
  r2_vals <- gdlite_diagnostics$r2_vals_by_k
  optimal_k <- gdlite_diagnostics$optimal_k
  png_file <- file.path(output_dir, "gdlite_r2_by_k.png")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  grDevices::png(png_file, width = 600, height = 400)
  k_vals <- seq_along(r2_vals)
  graphics::plot(k_vals, r2_vals, type = "b", col = "darkgreen", xlab = "Number of PCs (k)",
                 ylab = "Cross-validated R2", main = "GLMdenoise-Lite: R2 by k")
  if (!is.null(optimal_k) && optimal_k > 0 && optimal_k <= length(r2_vals)) {
    graphics::points(optimal_k, r2_vals[optimal_k], pch = 19, col = "red", cex = 1.5)
    graphics::text(optimal_k, r2_vals[optimal_k], labels = sprintf("k=%d", optimal_k), pos = 3, col = "red")
  }
  grDevices::dev.off()
  png_file
}
