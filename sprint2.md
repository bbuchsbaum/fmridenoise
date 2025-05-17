**Sprint Goal for Sprint 2: Implement Iterative ND-X Core, Introduce Auto-Adaptive Hyperparameters, and Develop "Annihilation Mode" & Enhanced Diagnostics.**
This sprint transitions ND-X from a single-pass system to its full iterative, self-tuning potential, directly targeting GLMdenoise's known limitations.

** always check off the tickets that are done. **

---
## ND-X Development: Sprint 2 Plan
**Sprint Duration:** (Estimate 3-4 weeks, given the complexity)
**Key Objectives:**
1. Implement the main iterative refinement loop (Pass 1, 2...).
2. Integrate full `NDx_FusedFIR` with k-medoids clustering and initial auto-adaptive cluster count logic.
3. Implement auto-adaptive rank for RPCA (`k_elbow`).
4. Implement BIC/AIC-based selection for spectral sines.
5. Implement the full Anisotropic Ridge regression with GCV-tuned lambdas for different subspaces.
6. Develop the "Annihilation Mode" by integrating GD-Lite PCs and orthogonalizing ND-X noise components.
7. Begin development of the HTML diagnostic report and JSON "Certified Clean" sidecar.
8. Introduce handling of the sparse component (`S`) from RPCA for spike TR down-weighting.
---
**Granular Tickets for Sprint 2:**
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
8. **NDX-8: [Workflow] Assemble Single-Pass ND-X Workflow** [X] 
* Create a main function `ndx_run_sprint1(Y_fmri, events, ...)` that calls NDX-2 through NDX-7 in sequence for a *single pass*.
* The `Y_residuals_current` from Pass-0 feeds into RPCA and Spectral.
* The HRFs from NDX-3 are used to build the task design for AR(2) fitting and Ridge.
* Nuisance regressors from RPCA and Spectral are included in the design for AR(2) and Ridge.
9. **NDX-9: [Diagnostics] Implement Initial Denoising Efficacy Score (DES)** [X]
* Implement `Calculate_DES(current_residuals_unwhitened, VAR_BASELINE_FOR_DES)`.
* Implement `Apply_Inverse_AR_Filter_To_Residuals`.
* Calculate and report DES after the single ridge regression pass.
10. **NDX-10: [Testing] Unit Tests & Basic Integration Test** [X] (Partially done for NDX-2)
* Write basic unit tests for key helper functions (e.g., FIR design builder, AR coefficient estimation).
* Create a small synthetic dataset and run the `ndx_run_sprint1` workflow to ensure parts connect and produce outputs of expected dimensions without errors.
**Epic: NDX-IterativeRefinement - Enabling Full Iterative Power & Auto-Adaptation**
[ ] 1. **NDX-11: [Workflow] Implement Main Iterative Loop & Convergence Logic**
* Refactor `ndx_run_sprint1` into `NDX_Process_Subject` (as per full pseudocode).
* Implement the `FOR pass_num FROM 1 TO MAX_PASSES` loop.
* Pass `Y_residuals_current` (unwhitened from previous pass) to RPCA and Spectral modules.
* Implement convergence checks based on `MIN_DES_GAIN_CONVERGENCE` and `MIN_RHO_NOISE_PROJECTION_CONVERGENCE`.
* Store `BETA_HISTORY_PER_PASS` and `DIAGNOSTICS_PER_PASS`.
[ ] 2. **NDX-12: [HRF-Est] Integrate Full `NDx_FusedFIR` with K-Medoids Clustering**
* Implement k-medoids clustering (e.g., using `cluster::pam`) on "good voxels" based on residuals from the *previous* ND-X pass (or Pass-0 for the first iteration).
* Fully integrate the `NDx_FusedFIR` R function (from the implementation-ready pseudocode), which includes robust FIR estimation.
* **Incorporate robustness for sparse events within FIR estimation logic:**
    * Implement event count check per condition (after spike masking).
    * Add user option `opts_hrf$hrf_min_events_for_fir` (e.g., default 6). If actual events < this, return a default zero/damped HRF or use a canonical HRF (strategy to be refined, simplest is to return zero/damped FIR).
    * Add user option `opts_hrf$hrf_low_event_threshold` (e.g., default 12) and `opts_hrf$hrf_target_event_count_for_lambda_scaling` (e.g., default 20).
    * If actual events are between `hrf_min_events_for_fir` and `hrf_low_event_threshold`, scale fused-lasso lambdas by `sqrt(hrf_target_event_count_for_lambda_scaling / actual_events)` before the final fit after CV.
* Implement initial `Auto_Adapt_HRF_Clusters` logic (e.g., merge highly similar, split high-variance).
* Ensure `motion_outlier_mask` (derived from `S` component of RPCA or FD, i.e., `current_spike_TR_mask`) and `voxel_R2` (from previous pass) are correctly passed and used for robust `ybar` calculation for clustering and per-cluster HRF estimation.
[X] 3. **NDX-13: [RPCA] Implement Auto-Adaptive Rank for RPCA-L**
* Implement `Auto_Adapt_RPCA_Rank` function using `k_elbow` logic (singular value drop-off, min/max clamping) as specified. (Function implemented and vetted as correct).
* Update the RPCA module (NDX-4) to use this adaptive rank in subsequent iterations. (Workflow calls adaptive function; input singular values for adaptation need to be correctly piped from previous pass RPCA results).
[X] 4. **NDX-14: [Spectral] Implement BIC/AIC Selection for Spectral Sines**
* Implement `Select_Significant_Spectral_Regressors` using `ΔBIC` (or AIC) to prune candidate spectral sines, as specified.
* Update the Spectral module (NDX-5) to use this selection logic.
[X] 5. **NDX-15: [Spikes-S] Implement Handling of Sparse Component `S` from RPCA**
* After RPCA, extract the `S` matrix. (Done in `ndx_rpca_temporal_components_multirun`)
* Implement `S_t_*` calculation (median absolute value across voxels per TR) to generate `spike_TR_mask` (global TR flag). (Implemented in `ndx_rpca_temporal_components_multirun`)
* Sprint 2 focus: Implement `spike_TR_mask` and ensure it's used in `NDx_FusedFIR`'s `w_TR`. (Mask generated and pipeline uses it for HRF estimation via `current_spike_TR_mask`)
* Start exploring how `S_t,v` can influence precision for the ridge. (Exploratory)
* (Note: Include an early abort if `sum(abs(S))==0` to skip mask generation.) (Implemented per run)
* (Follow-up: `ndx_rpca_temporal_components_multirun` should also return `V_global_singular_values` for NDX-13 and optionally `S_matrix_cat` for diagnostics/future use.)
[ ] 6. **NDX-16: [AnisoRidge] Implement Full Anisotropic Ridge Regression**
* Refactor the basic ridge (NDX-7) to `ndx_solve_anisotropic_ridge`.
* Implement creation of projection matrices: `P_GD`, `P_Unique` (if Annihilation Mode), `P_Noise` (if not Annihilation Mode), and `P_Signal`.
* Implement GCV for tuning `lambda_parallel` (and `lambda_unique_nuisance` if applicable) and `lambda_signal_friendly_perp`. Ensure lambdas are scaled by `res_var_whitened_estimate`.
* Implement `Update_Lambda_Aggressiveness` based on `rho_noise_projection`. (Note: Log `lambda_parallel` (λ∥) per pass to `DIAGNOSTICS_PER_PASS`.)
[ ] 7. **NDX-17: [Annihilation] Implement "Annihilation Mode" Core Logic**
* Integrate `U_GD_PCs` (from a quick GD-Lite run, NDX-Ext-1 if made external, or internal call).
* Implement `Gram_Schmidt_Orthogonalize_Against(U_NDX_Nuisance_Combined, U_GD_PCs)` to get `U_NDX_Unique_Nuisance`.
* Ensure Anisotropic Ridge (NDX-16) can handle the three-projector (`P_GD`, `P_Unique`, `P_Signal`) scenario.
**Epic: NDX-Reporting - Initial Diagnostics & Certification**
[ ] 8. **NDX-18: [Diagnostics] Develop Initial HTML Diagnostic Report Structure**
* Design the basic layout for the HTML report.
* Implement functions to generate and embed initial plots:
    * DES per pass.
    * PSD of residuals (Pass-0 vs. final ND-X pass).
    * Spike carpet plot from `S` matrix (if `S` handling is sufficiently mature).
    * Display key adaptive hyperparameter choices.
[ ] 9. **NDX-19: [Certificate] Implement "ND-X Certified Clean" JSON Sidecar Generation**
* Define the JSON structure.
* Populate with initial fields: `ndx_version`, final `DES`, `num_passes_converged`, final adaptive hyperparameter values, final `rho_noise_projection`.
* Write function to save this JSON file.
[ ] 10. **NDX-20: [Testing] Iterative Workflow & Annihilation Mode Integration Tests**
* Expand tests from Sprint 1 to cover the iterative workflow with 2-3 passes.
* Test the Annihilation Mode path specifically: ensure `U_GD_PCs` are generated, orthogonalization occurs, and the ridge uses the correct projectors.
* Verify convergence logic and adaptive hyperparameter updates.
**External Dependencies/Assumptions for Sprint 2:**
* `glmgen::fusedlasso` is available and its parameter mapping is understood.
* RPCA R package (`rpca` or similar) is functional.
* Multitaper spectrum package (`multitaper` or `psd`) is functional.
* Clustering package (`cluster::pam`) is functional.

**Out of Scope for Sprint 2 (deferred to Sprint 3):**
* Full, polished HTML report with all progressive enhancement visualizations and Annihilation Verdict grading.
* Final, highly optimized `n_threads` usage across all C++/Rcpp components.
* Exhaustive stress-testing against diverse "pathological" datasets (initial robustness via fallbacks in place).
* Advanced `S_t,v` based precision re-weighting fully integrated into the ridge objective.
* User-configurable manual overrides for all auto-adaptive parameters (focus on auto-pilot first).
* Hierarchical (partial-pool) shrinkage of HRF estimates towards a grand average (enhancement for `NDx_FusedFIR`).

---
## Addendum: Detailed Plan for NDX-15 (Spike Handling from RPCA 'S' Component)

This addendum outlines the specific implementation steps for generating and utilizing the `spike_TR_mask` as part of ticket NDX-15.

**A. Modifications to `ndx_rpca_temporal_components_multirun` (in `R/ndx_rpca.R`):**

1.  **Initialization (before per-run loop):**
    *   Initialize `per_run_spike_TR_masks <- list()` to store logical masks for each run.

2.  **Inside the Per-Run RPCA Loop (for each `run_name` and its sparse component `S_r_t` which is Voxels x Time_r):
    *   **a. Early Abort for Empty S:** If `sum(abs(S_r_t), na.rm = TRUE) < 1e-9`, then `spike_TR_mask_run <- rep(FALSE, ncol(S_r_t))`. Store this in `per_run_spike_TR_masks[[run_name]]` and proceed to the next run.
    *   **b. Calculate Per-TR Spike Summary (`s_t_star_run`):**
        *   `s_t_star_run <- apply(abs(S_r_t), 2, stats::median, na.rm = TRUE)`.
    *   **c. Threshold `s_t_star_run` to get `spike_TR_mask_run`:**
        *   Retrieve/define `rpca_spike_mad_thresh` and `rpca_spike_percentile_thresh` from `current_opts` (add to defaults in `ndx_rpca_temporal_components_multirun`).
        *   Calculate `mad_s_t_star`.
        *   If `mad_s_t_star` is too small, use percentile threshold; otherwise, use MAD-based threshold.
        *   `spike_TR_mask_run <- s_t_star_run > threshold_s_t_star`.
    *   **d. Store Per-Run Mask:** `per_run_spike_TR_masks[[run_name]] <- spike_TR_mask_run`.

3.  **After the Per-Run RPCA Loop (before returning from `ndx_rpca_temporal_components_multirun`):
    *   **e. Combine Per-Run Masks into Global `spike_TR_mask`:**
        *   `global_spike_TR_mask <- unlist(per_run_spike_TR_masks[names(Y_residuals_list)], use.names = FALSE)`.
        *   Validate length against `nrow(Y_residuals_cat)`.
    *   **f. Update Function Return Value:** Ensure the function returns a list including `spike_TR_mask = global_spike_TR_mask` (along with `C_components`, `S_matrix_cat` if needed, and `V_global_singular_values`).

**B. Modifications to `NDX_Process_Subject` (in `R/ndx_workflow.R`):**

1.  **Initialize `current_spike_TR_mask`:** Existing initialization is good.
2.  **Process `rpca_out`:** Ensure `rpca_out$spike_TR_mask` is correctly extracted and ORed with `current_spike_TR_mask`.
3.  **Pass `current_spike_TR_mask` to `ndx_estimate_initial_hrfs`:** Already done.

**C. Add New User Options (to `default_opts` in `ndx_rpca_temporal_components_multirun`):**

*   `rpca_spike_mad_thresh` (e.g., default `3.0`).
*   `rpca_spike_percentile_thresh` (e.g., default `0.98`).

**D. Testing (as part of NDX-20):**

*   Test cases with known spikes.
*   Verify `spike_TR_mask` generation and its effect in HRF estimation.