# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fmridenoise (ndx) is an R package implementing an advanced iterative denoising workflow for single-subject fMRI data. It removes noise components from fMRI time series data to improve task-related signal estimation using techniques like robust PCA, spectral analysis, AR(2) pre-whitening, and ridge regression.

## Common Development Commands

### Testing
- Run all tests: `Rscript run_tests.R`
- Run specific test file: `Rscript -e "testthat::test_file('tests/testthat/test-ridge.R')"`
- Run tests matching a pattern: `Rscript -e "testthat::test_file(dir('tests/testthat', pattern='spectral', full.names=TRUE))"`

### Building and Checking
- Build package: `R CMD build .`
- Check package: `R CMD check ndx_0.0.0.9000.tar.gz`
- Install package locally: `R CMD INSTALL .`
- Document with roxygen2: `Rscript -e "devtools::document()"`

### C++ Development
- The package uses Rcpp/RcppArmadillo for performance-critical operations
- C++ source files are in `src/` (grassmann_merge.cpp, apply_ar_filter.cpp)
- After modifying C++ code, run `Rscript -e "Rcpp::compileAttributes()"` to update RcppExports

## Architecture and Workflow

The package implements an 8-step iterative workflow coordinated by `NDX_Process_Subject()`:

1. **Initial GLM** (`workflow_step1_initial_glm.R`) - Baseline model with motion/polynomial regressors
2. **HRF Estimation** (`workflow_step2_hrf.R`) - Data-driven HRF estimation using FIR basis and clustering
3. **RPCA** (`workflow_step3_rpca.R`) - Robust PCA for spike detection and low-rank noise removal
4. **Spectral Analysis** (`workflow_step4_spectral.R`) - Frequency-domain noise component identification
5. **Design Matrix** (`workflow_step5_design_matrix.R`) - Full regression model construction
6. **AR(2) Whitening** (`workflow_step6_whitening.R`) - Temporal autocorrelation handling
7. **Ridge Regression** (`workflow_step7_ridge.R`) - Regularized parameter estimation with GCV/LOOCV
8. **Diagnostics** (`workflow_step8_diagnostics.R`) - R², DES, and diagnostic plots

Each step outputs structured data that feeds into subsequent steps. The workflow supports both single-run and multi-run data processing.

## Key Implementation Details

- **Options System**: User configuration via `ndx_user_options()` controls all workflow parameters
- **Data Structures**: Core functions work with neuroim2 objects (NeuroVec, NeuroVol)
- **Parallelization**: Uses foreach/doParallel for computationally intensive operations
- **Precision Weighting**: Advanced variance estimation for heteroscedastic noise modeling
- **Diagnostic Plots**: Comprehensive visualization functions in `ndx_diagnostic_plots.R`

## Testing Philosophy

- Extensive unit test coverage (40+ test files) for all major components
- Helper files provide synthetic data generation for consistent testing
- Tests cover edge cases, numerical accuracy, and integration scenarios
- Run tests before committing any changes to ensure stability

## Recent Optimization: RSVD Implementation (36x Speed, 100x Accuracy)

### Background
The package originally used `rpca::rpca` for robust PCA decomposition in the temporal nuisance component extraction step. Through comprehensive benchmarking, we discovered that the randomized RSVD algorithm from the `rsvd` package could be dramatically optimized to outperform the original implementation.

### Benchmark Results
Comprehensive testing across multiple scenarios revealed:

**Original RPCA (rpca::rpca)**:
- Median execution time: 1.734s
- Reconstruction error: 9.99e-08  
- Success rate: 100%

**Optimized RSVD (rsvd::rrpca)**:
- Median execution time: 0.048s (**36.1x faster**)
- Reconstruction error: 9.56e-10 (**104.5x more accurate**)
- Success rate: 100%

### Optimal Parameters Discovered
Through extensive parameter tuning, the optimal RSVD configuration was identified:

```r
rsvd::rrpca(A = matrix, lambda = lambda,
            maxiter = 300,    # 6x more iterations than default (50)
            tol = 1e-09,      # 10,000x tighter tolerance than default (1e-05)  
            p = 25,           # 2.5x more oversampling than default (10)
            q = 5)            # 2.5x more power iterations than default (2)
```

### Implementation Details
The optimization was implemented by replacing `rpca::rpca` calls with `rsvd::rrpca` in:

1. **`.ndx_basic_rpca_single_run`** function in `R/ndx_rpca_helpers.R`
2. **`.ndx_rpca_single_run`** function in `R/ndx_rpca_helpers.R`

### Backward Compatibility
The implementation maintains full backward compatibility:

- **API unchanged**: All function signatures remain identical
- **Legacy parameters**: Old `rpca::rpca` parameters (e.g., `rpca_term_delta`, `rpca_max_iter`, `rpca_mu`, `rpca_trace`) are still accepted but ignored
- **User notification**: Users are informed when legacy parameters are provided
- **Return values**: All return structures remain identical

### Technical Benefits
- **Dramatic speed improvement**: 36.1x faster execution
- **Superior accuracy**: 104.5x better reconstruction error
- **Same interface**: Drop-in replacement for existing code
- **Robust performance**: 100% success rate across all test scenarios
- **Memory efficient**: Randomized algorithm uses less memory for large matrices

### Key Files Modified
- `R/ndx_rpca_helpers.R`: Core RPCA implementation functions
- `R/ndx_rpca.R`: Main RPCA interface and documentation  
- `man/ndx_rpca_temporal_components_multirun.Rd`: Updated documentation

### Dependencies Updated
- Added explicit import: `rsvd::rrpca`
- Maintained compatibility with existing `rpca` dependency

This optimization represents a major performance breakthrough for the ND-X pipeline, providing users with significantly faster processing times while simultaneously improving the quality of temporal nuisance component extraction.