Given that ND-X will be implemented in R and aims to "beat the living tar out of GLMdenoise," and we'll leverage `fmrireg` for HRF convolution modeling, let's define Sprint 1 for ND-X development with granular tickets.
**Sprint Goal for Sprint 1: Establish Core ND-X Pipeline with Initial Denoising Layers & Foundational HRF Modeling (using `fmrireg`).**
This sprint focuses on getting the basic structure of ND-X in place, implementing the first pass of several key modules, and ensuring robust data handling and initial HRF estimation. The auto-adaptive features and full annihilation mode will be built upon this foundation in later sprints.

root folder is  "fmridenoise" not "ndx"
---
## ND-X Development: Sprint 1 Plan
**Sprint Duration:** (Estimate 2-3 weeks, depending on team size/focus)
**Key Objectives:**
1. Implement "Pass 0" for initial residual generation. [X]
2. Integrate `fmrireg` for flexible, condition-specific FIR HRF estimation (non-clustered, basic fused-lasso for now). [X]
3. Implement the first pass of Robust PCA (RPCA-L) for low-rank noise component identification. [X]
4. Implement the first pass of Multi-taper Spectral analysis for periodic nuisance regressor identification. [ ]
5. Develop the initial AR(2) pre-whitening module. [ ]
6. Implement a basic (isotropic or simplified anisotropic) ridge regression. [ ]
7. Set up basic data loading, preprocessing, and output structures. [X] (Covered by NDX-1 utils)
8. Establish initial diagnostic metrics (e.g., basic DES). [ ]
---
**Granular Tickets for Sprint 1:**
**Epic: NDX-CorePipeline - Establish Foundational ND-X Workflow**
1. **NDX-1: [Setup] Project Setup & `fmrireg` Integration** [X]
* Create R package structure for ND-X. [X]
* Define core data structures/classes for managing fMRI data, design matrices, and ND-X parameters. (Placeholder created) [X]
* Ensure `fmrireg` is a dependency and basic `fmrireg` functions for HRF object creation and design matrix generation can be called. [X]
* Write utility functions for loading fMRI data (e.g., NIfTI via `neuroim` or `oro.nifti`) and basic pre-processing (demeaning, detrending if not part of Pass 0). (Placeholder created) [X]
2. **NDX-2: [Pass-0] Implement Initial OLS Residual Generation** [X]
* Implement function `ndx_pass0_residuals(Y_fmri, events, motion_params, run_idx, TR)`. [X]
* Use `fmrireg::hrf_spmg1()` or similar for canonical HRF. [X]
* Use `fmrireg` utilities (e.g., `event_factor`, `covariate`, `baseline`) to build the `X_basic_task` and `X_basic_nuisance` design matrices. [X]
* Handle run dummies if `Detect_Run_Specific_Tasks(events)` indicates task differences. (Deferred - not explicitly implemented yet, current model is simpler)
* Perform OLS and return `Y_residuals_current` and `VAR_BASELINE_FOR_DES`. [X]
3. **NDX-3: [HRF-Est] Initial FIR HRF Estimation with `fmrireg` & Fused Lasso** [X]
* Implement `ndx_estimate_initial_hrfs(Y_fmri, events, TR, R2_pass0, spike_TR_mask, user_options)`. [X]
* For Sprint 1, simplify:
    * No clustering yet; estimate one "global" or "good-voxel-average" HRF per condition. [X]
    * Use `fmrireg::build_FIR_design` (or the one provided in `NDx_FusedFIR` pseudocode - used `fmrireg::event_model` with tent basis). [X]
    * Use `glmgen::fusedlasso` directly on `ybar` (robust mean from good voxels) for each condition's `X_q`. [X]
    * Implement basic GCV for `lambda`, `gamma` for these initial HRFs. [X]
    * Implement `project_cone` for physiological plausibility. (Basic version implemented) [X]
* Output: A tibble of HRF impulse responses. [X]
* *Note: Full `NDx_FusedFIR` with clustering will be in a later sprint.*
4. **NDX-4: [RPCA] Implement Robust PCA (L-component) on Residuals** [X]
* Implement a function (e.g., `ndx_rpca_temporal_components_multirun`) that takes `Y_residuals` (concatenated) and `run_idx`. [X]
* This function will implement the per-run RPCA strategy detailed in the `proposal.md` Addendum:
    * For each run, apply RPCA to `t(Y_residuals_run)` (voxels x time) to get voxel-space components `V_r`. [X]
    * Merge these `V_r` components to obtain `V_global`. For Sprint 1, the primary merge strategy will be "concatenate and SVD" (`V_all_concat = do.call(cbind, V_list_valid); V_global = svd(V_all_concat)$u`), as it's mathematically equivalent to the Grassmann mean under equal weighting and often faster if memory permits. An option for iterative Grassmann averaging will also be included for flexibility and scalability. [X]
    * Project per-run residuals onto `V_global` to get run-specific temporal nuisance regressors `C_r`. [X]
* `k_global_target` will define the number of components in `V_global`. [X]
* The function should return the concatenated `C_r` matrix. [X]
* Implement basic "glitch ratio" check (can be per-run and summarized) and rank capping safeguards. [X]
* *Decision: Implement per-run RPCA. For merging `V_r` components into `V_global`, Sprint 1 will offer both the "concatenate and SVD" method (default, fast for moderate data, valid Grassmann mean) and an iterative Grassmann averaging method (memory-efficient fallback). This aligns with the `proposal.md` addendum discussion on RPCA strategies and their trade-offs.* [X]
* Write unit tests for `ndx_rpca_temporal_components_multirun`. [X]
5. **NDX-5: [Spectral] Implement Multi-Taper Spectral Nuisance Identification** [X]
* Implement `ndx_spectral_sines(mean_residual_for_spectrum, TR, n_sine_candidates, nyquist_guard_factor)`. [X]
* Use R's `spec.mtm` (from `multitaper` package) or `psd::pspectrum` for multi-taper spectrum estimation. [X]
* Identify peaks and generate sine/cosine regressors `U_Spectral_Sines`. [X]
* For Sprint 1, selection can be based on peak prominence; full BIC/AIC filtering later. [X]
6. **NDX-6: [Whitening] Implement AR(2) Pre-whitening** [X]
* Implement `ndx_ar2_whitening(Y_data, X_design_full, Y_residuals_for_AR_fit)`. [X]
* Estimate AR(2) coefficients voxel-wise from `Y_residuals_for_AR_fit` (e.g., using `stats::ar.yw`). [X]
* Implement functions `Apply_AR_Filter_To_Data` and `Apply_AR_Filter_To_Design`. [X]
* Return `Y_whitened`, `X_whitened`, and `AR_coeffs`. [X]
7. **NDX-7: [Ridge] Implement Basic Ridge Regression** [X]
* Implement `ndx_solve_ridge(Y_whitened, X_whitened, lambda_ridge)`. [X]
* For Sprint 1, this can be a standard isotropic ridge. `lambda_ridge` can be selected via a simple GCV on a small grid. [X]
* Return `betas_whitened`. [X]
* Implement `Extract_And_Unwhiten_Task_Betas` (initial version, may need refinement later). [X]
8. **NDX-8: [Workflow] Assemble Single-Pass ND-X Workflow** [ ]
* Create a main function `ndx_run_sprint1(Y_fmri, events, ...)` that calls NDX-2 through NDX-7 in sequence for a *single pass*.
* The `Y_residuals_current` from Pass-0 feeds into RPCA and Spectral.
* The HRFs from NDX-3 are used to build the task design for AR(2) fitting and Ridge.
* Nuisance regressors from RPCA and Spectral are included in the design for AR(2) and Ridge.
9. **NDX-9: [Diagnostics] Implement Initial Denoising Efficacy Score (DES)** [ ]
* Implement `Calculate_DES(current_residuals_unwhitened, VAR_BASELINE_FOR_DES)`.
* Implement `Apply_Inverse_AR_Filter_To_Residuals`.
* Calculate and report DES after the single ridge regression pass.
10. **NDX-10: [Testing] Unit Tests & Basic Integration Test** [ ] (Partially done for NDX-2)
* Write basic unit tests for key helper functions (e.g., FIR design builder, AR coefficient estimation).
* Create a small synthetic dataset and run the `ndx_run_sprint1` workflow to ensure parts connect and produce outputs of expected dimensions without errors.
**Out of Scope for Sprint 1 (deferred to later sprints):**
* Full iterative refinement loop (Pass 1, 2...).
* Advanced auto-adaptive hyperparameter tuning (k_RPCA elbow, HRF cluster merging/splitting, BIC for sines, adaptive lambda aggressiveness).
* Full `NDx_FusedFIR` with k-medoids clustering.
* Full Anisotropic Ridge regression with multiple projectors and GCV for `λ_parallel`, `λ_perp`.
* Annihilation Mode (GD-Lite, orthogonalization, three-block ridge).
* Comprehensive HTML diagnostic report and JSON "Certified Clean" sidecar.
* Spike (`S` component from RPCA) handling beyond basic outlier masking for HRF estimation.
* Advanced parallelization (`n_threads` propagation to all Rcpp/parallelizable steps).
---
This Sprint 1 plan focuses on building the critical scaffolding and first versions of each core module. It leverages `fmrireg` for its strengths in HRF/design handling and `glmgen` for a specific solver, allowing the ND-X team to focus on the novel denoising logic. Subsequent sprints will layer on the iterative refinement, auto-adaptation, advanced regularization, and the "annihilation" features.