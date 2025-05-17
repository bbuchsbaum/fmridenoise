**Sprint Goal for Sprint 3: Achieve Production-Ready ND-X: Full Diagnostics, Robustness Hardening, Performance Optimization, and Definitive Validation.**
This sprint is about ensuring ND-X is not just powerful, but also reliable, user-friendly, transparent, and its superiority is irrefutably demonstrated.
---
## ND-X Development: Sprint 3 Plan
**Sprint Duration:** (Estimate 3-4 weeks, focused on polish, testing, and documentation)
**Key Objectives:**
1. Finalize and polish the comprehensive HTML diagnostic report, including progressive enhancement visualizations and the graded "Annihilation Verdict."
2. Complete the "ND-X Certified Clean" JSON sidecar with all relevant metrics.
3. Implement advanced, robust handling of the sparse component `S` from RPCA, fully integrating `S_t,v`-based precision re-weighting.
4. Conduct thorough stress-testing across diverse and "pathological" datasets, refining auto-adaptive fallbacks and guardrails.
5. Optimize performance, including systematic `n_threads` utilization and profiling critical code paths.
6. Perform definitive benchmark comparisons against GLMdenoise on multiple public datasets.
7. Finalize user documentation, tutorials, and packaging for release.
---
**Granular Tickets for Sprint 3:**
**Epic: NDX-PolishAndValidate - Finalizing for Release and Demonstrating Dominance**
1. **NDX-21: [Diagnostics] Finalize HTML Diagnostic Report with Progressive Enhancements**
* Implement the four-panel slider visualization showing DES, β-stability, and residual PSD at stages: (1) GD-equivalent, (2) + ND-X unique regressors, (3) + anisotropic ridge, (4) + iterative refinement (final).
* Include the graded "Annihilation Verdict" (Tie, Win, Decisive Win, Annihilation) based on `Var_NDX_unique / Var_GD_only` ratio and/or percentage of total removable nuisance.
* Ensure all key plots (PSD ratios, spike carpet, beta stability, Ljung-Box for whiteness) are polished and clearly presented.
* Add a section summarizing final adaptive hyperparameter choices and convergence information.
2. **NDX-22: [Certificate] Complete "ND-X Certified Clean" JSON Sidecar**
* Ensure all specified fields are present and accurately populated: `ndx_version`, `DES`, `var_ratio` (for Annihilation Verdict), `ljung_box_p`, final adaptive hyperparams (`k_rpca`, `n_hrf_clusters`, `n_spectral_regressors`), `passes_converged`, `verdict`.
* Add a timestamp and potentially a hash of key input parameters for provenance.
3. **NDX-23: [Spikes-S] Advanced `S_t,v`-Based Precision Reweighting**
* Fully integrate `S_t,v` (sparse component values at specific voxel-time points) into the precision estimation.
* This means that the voxel-wise residual variance used for the AR(2) model and, crucially, for scaling the anisotropic ridge penalties (and potentially as weights in the ridge objective itself if using a weighted least squares formulation) is directly informed by `S_t,v`.
* If `S_{t,v}` is large, the precision `w_{t,v}` should be very low for that specific point, effectively down-weighting its influence in the GLM fit.
* Ensure this interacts correctly with the AR(2) whitening process.
4. **NDX-24: [Robustness] Exhaustive Stress-Testing & Fallback Refinement**
* Test ND-X on a diverse set of fMRI datasets, including:
* High-motion subjects.
* Data with severe scanner artifacts (e.g., known spike issues).
* Low SNR / subtle task activation datasets.
* Datasets with very few runs or short runs.
* Datasets with unusual TRs or event timings.
* Monitor the behavior of auto-adaptive rules and guardrails (NDX-Ext-Pathology from prior discussions).
* Refine fallback mechanisms (e.g., `rpca-max-rank`, `hrf-fallback-canonical`, `force-spectral N` options and their auto-triggers) based on these tests to ensure ND-X *never* performs worse than a basic GLM or GLMdenoise.
5. **NDX-25: [Optimization] Performance Profiling & Parallelization**
* Profile the entire ND-X pipeline on representative datasets to identify computational bottlenecks.
* Systematically implement and test parallelization using `n_threads` (e.g., via `RcppParallel` for custom C++ code, or ensuring R package calls leverage available cores) for major computational blocks (RPCA, HRF estimation loops, GCV for ridge, etc.).
* Optimize R code for memory efficiency and speed (vectorization, avoiding unnecessary copies).
* Aim to consistently meet or beat the "~9 min on 12-cores" target for typical datasets.
6. **NDX-26: [Benchmarking] Definitive GLMdenoise Comparison**
* Select 2-3 diverse, public fMRI datasets (e.g., HCP 7T, NSD, a clinical dataset if available).
* Run both ND-X (auto-pilot, Annihilation Mode with `--benchmark-gd` flag) and a fully tuned, official GLMdenoise implementation on these datasets.
* Collect and compare all key metrics: DES, Annihilation Verdict scores (`Var_NDX_unique`, `Var_GD_only`), cross-run RSA, leave-run-out decoding accuracy, group-level statistical power (if feasible with simulated or appropriate real group data), and wall-time.
* Prepare tables and figures summarizing these results for publication/documentation.
**Epic: NDX-ReleasePrep - Documentation & Packaging**
7. **NDX-27: [Docs] Comprehensive User Documentation & Tutorials**
* Write detailed documentation for installing and running ND-X.
* Explain all user-configurable options and their defaults.
* Provide guidance on interpreting the HTML diagnostic report and the "Certified Clean" JSON.
* Create tutorials/vignettes with example datasets, including one demonstrating the "Progressive Enhancement Visualization."
* Document the underlying algorithms and key citations.
8. **NDX-28: [Packaging] Prepare R Package for Distribution**
* Ensure the R package passes `R CMD check` cleanly.
* Include all necessary dependencies and version requirements.
* Prepare for potential submission to CRAN or distribution via GitHub.
* Consider creating a Docker image for easy, reproducible deployment.
9. **NDX-29: [API] Finalize API for Downstream Integration (e.g., CV-PWS)**
* Confirm the structure of output betas, precision maps, and the JSON sidecar is stable and convenient for tools like CV-PWS to consume.
* Provide example code snippets for how CV-PWS (or other tools) could read the "Certified Clean" certificate and potentially adapt its behavior.
10. **NDX-30: [Showcase] Prepare Pre-print / Publication Materials**
* Draft a manuscript detailing ND-X's methodology, features, and benchmark results.
* Create figures and tables based on NDX-26 and the diagnostic report features.
* Highlight the "Annihilation Mode" and "Certified Clean" concepts as key innovations.