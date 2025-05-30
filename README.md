# ND-X: Next-Generation fMRI Denoising

ND-X implements an iterative denoising workflow for single-subject fMRI data. The package integrates robust PCA, spectral analysis, AR(2) pre-whitening and ridge regression to automatically identify nuisance components and refine task estimates. It aims to surpass existing methods such as GLMdenoise by combining multiple data-adaptive steps.

## Key Functions

- `NDX_Process_Subject()` – run the full denoising workflow for one subject.
- `ndx_default_user_options()` – return a list of default workflow options.
- `Auto_Adapt_RPCA_Rank()` – choose an appropriate RPCA rank from singular values.
- `Select_Significant_Spectral_Regressors()` – pick sine/cosine regressors based on information criteria.
- `calculate_DES()` – compute the Denoising Efficacy Score.
- `ndx_generate_html_report()` – create an HTML report of diagnostics.

See individual function documentation for many other helpers (HRF estimation, whitening, ridge solving, etc.).

## Basic Usage

Inputs to `NDX_Process_Subject()` are matrices of fMRI data and motion regressors, an event table and run index, the TR and an optional spike mask. Most behaviour is controlled through the `user_options` list; `ndx_default_user_options()` gives a sensible starting point.

### Quick Start Example

```r
library(ndx)

# fmri_data, events_df, motion_matrix and run_index should be prepared by the user
opts <- ndx_default_user_options()
results <- NDX_Process_Subject(
  Y_fmri       = fmri_data,
  events       = events_df,
  motion_params = motion_matrix,
  run_idx      = run_index,
  TR           = 2,
  user_options = opts
)
```

The returned list contains estimated HRFs, nuisance components, final betas and diagnostic metrics. Diagnostic plots can be written with `ndx_generate_html_report(results, results$pass0_residuals, TR = 2)`.

## Testing

Unit tests for the package are located under `tests/` and can be run with

```r
Rscript run_tests.R
```

## License

This project is released under the MIT license.
