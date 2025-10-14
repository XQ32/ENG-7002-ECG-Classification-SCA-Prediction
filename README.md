## Project Overview

This repository implements a complete single-lead ECG pipeline on SHHS1:

- Heartbeat detection and beat-type recognition (PVC vs Other)
- Record-wise generation of PVC predictions and derived info (`*_info.mat`) without human labels
- A variable-window feature system of “PVC → sinus recovery dynamics” to build an SCA risk model (binary: Alive=1 / Dead=0)

Key scripts:

- `generate_beats_classifier.m`: extract per-beat features from rpoints-annotated records and train the beat-type (PVC/Other) classifier
- `predict_shhs1_all.m`: iterate over all EDFs and produce `*_info.mat` for each record using the trained beat-type model
- `generate_sca_classifier_v3.m`: read `*_info.mat`, compute record-level features from PVC-triggered recovery, and train/evaluate the SCA model (v3)

Support layers:

- Preprocessing: `ecgFilter.m`
- Beat detection/labeling: `detectAndClassifyHeartbeats.m` (includes QRS helpers, T-wave seeds, and SQI gating via `assessBeatsQuality.m`)
- Per-beat features: `extractHeartbeatFeatures.m`
- `mcode/`: WFDB toolbox (PhysioNet) for MATLAB and examples

Outputs (v3):

- `results/post_ectopic_features_v3.mat`: original record-level feature table (with NaN)
- `results/SCA_trainingFeatureTable_v3.mat`, `results/SCA_testingFeatureTable_v3.mat`
- `results/sca_classifier_v3.mat`: model package (with importances and univariate effectiveness)
- `results/*.png`: ROC/PR/calibration/decision curves, etc.
- `results/*_v3.csv`: univariate effectiveness, ablation, SHAP importance

Recommended MATLAB version: R2021a+ (to support the SHAP computation path). Toolboxes required: Statistics and Machine Learning, Signal Processing.

Data layout requirement (relative to project root):

- `shhs/polysomnography/edfs/shhs1/*.edf`
- `shhs/polysomnography/annotations-rpoints/shhs1/*-rpoint.csv` (only when training the beat-type model)
- `shhs/datasets/shhs-cvd-summary-dataset-0.21.0.csv` (mapping `nsrrid → vital`)

---

## Quick Start

1) Train (or prepare) the beat-type (PVC/Other) model (optional)

- If `results/trainedClassifier_latest.mat` already exists, skip.
- Otherwise run: `generate_beats_classifier.m`
	- It will read rpoints-annotated records → filter → per-beat features → stratified split → train AdaBoost/tree ensembles → save `results/trainedClassifier_latest.mat`.

2) Generate `*_info.mat`

- Run: `predict_shhs1_all.m`
	- For each EDF: filtering → beat detection and features → use the model from step 1 to classify PVC → write required fields into `{record}_info.mat` in the same folder (see “latestRequiredFields” in the script).

3) Train the SCA model (v3)

- Run: `generate_sca_classifier_v3.m`
	- Iterate over `*_info.mat` → compute recovery-dynamics features for each PVC → aggregate to record level → (impute missing + MinMax normalization) → group-wise stratified split (by patient) → univariate evaluation → model training (AdaBoostM1 with threshold moving / Platt calibration) → evaluation and plots → save the model package.

Note: by default, Alive=1 / Dead=0 (`patientVital` is the response).

---

## Signal Processing and Detection (Brief)

- Preprocessing (`ecgFilter.m`):
	- Power-line notch (optional 50/60 Hz and 2nd/3rd harmonics); three baseline removal strategies (FIR / IIR / enhanced pipeline); low-pass; method 3 adds SG smoothing and Pan–Tompkins-style QRS enhancement.
- Beat detection (`detectAndClassifyHeartbeats.m`):
	- Absolute-amplitude + derivative-energy candidates → map back to |x| peaks; remove T-as-R errors; Q/S localization (energy-gated, polarity-adaptive, PVC-wide windows); QRS onset/offset (derivative-threshold tracking); optional template-correlation SQI gating (`assessBeatsQuality.m`).
- Per-beat features (`extractHeartbeatFeatures.m`):
	- RR before/after, R/Q/S amplitudes, QRS duration/area, with fallbacks and normalization.

These steps are orchestrated in `predict_shhs1_all.m` to produce the compact fields required by the v3 generator in `*_info.mat`.

---

## v3 Core: Full List and Explanation of “Features Used for Training”

`generate_sca_classifier_v3.m` computes many PVC-level recovery metrics and aggregates them to record-level features, then selects 15 core features based on experience and SHAP/univariate/AUC criteria. Below are the 15 features actually used for training (original name → renamed) and their meaning. We use dev to denote relative deviation: $\text{dev}_{RR}=|RR-\mu_{RR}|/\mu_{RR}$.

1) PVC_FreqPerHour (orig: PVCs_per_hour → new: PVC_FreqPerHour)
- Definition: PVC count per hour (normalized by total record duration).
- Meaning: PVC burden/frequency.

2) RRI_Recovery_TimeConst_Median (orig: tau_rr_median → new: RRI_Recovery_TimeConst_Median)
- Definition: the time constant $\tau$ from a log-linear fit on the first K points (K=8 by default) of dev_RR; aggregated by median over PVCs.
- Meaning: scale of exponential recovery speed; larger means slower recovery.

3) RRI_PostPVC_OscillAmp_Median (orig: oscill_rr_median → new: RRI_PostPVC_OscillAmp_Median)
- Definition: normalized oscillation index from sign flips of the first difference of dev_RR (`local_oscillation_index`); aggregated by median.
- Meaning: oscillation/underdamped behavior during recovery.

4) HRT_TurbOnset_AbnormalFrac (orig: HRT_TO_abnormal_frac → new: HRT_TurbOnset_AbnormalFrac)
- Definition: fraction of abnormal heart rate turbulence TO (threshold TO≥0 in the script).
- Meaning: abnormal vagal/sympathetic response to PVC (first post-PVC RR change).

5) RRI_PostPVC_OscillAmp_Log (orig: oscill_log → new: RRI_PostPVC_OscillAmp_Log)
- Definition: log transform of feature 3, log(max(oscill, eps)).
- Meaning: compress high oscillation values, reduce long-tail effects.

6) PVC_FreqPerHour_Sqrt (orig: pvc_freq_sqrt → new: PVC_FreqPerHour_Sqrt)
- Definition: square root transform of PVC frequency.
- Meaning: mitigate scale effects in very high-frequency records.

7) RRI_Oscill_x_PVC_Freq (orig: oscill_x_pvc_freq → new: RRI_Oscill_x_PVC_Freq)
- Definition: oscillation index × PVC frequency.
- Meaning: combined burden for high PVC frequency and strong oscillation.

8) RRI_Recovery_EarlyLate_Ratio (orig: early_late_ratio → new: RRI_Recovery_EarlyLate_Ratio)
- Definition: ratio of early to mid-term recovery half-life: `halflife_rr10_median / halflife_rr30_median`.
- Meaning: relational speed between early and mid recovery (>1 often indicates “fast early, slower mid”).

9) RRI_Recovery_TimeVariability_CV (orig: halflife_cv → new: RRI_Recovery_TimeVariability_CV)
- Definition: coefficient of variation (CV) of (50%) recovery half-life across PVCs.
- Meaning: instability/inconsistency of recovery times across PVCs.

10) Clinical_CompositeRisk_Score (orig: composite_risk_score → new: Clinical_CompositeRisk_Score)
- Definition: experience-weighted composite: oscill_norm(0.40) + pvc_norm(0.30) + hl10_norm(0.30), where each is normalized to [0,1] by common thresholds.
- Meaning: clinically interpretable combined risk indicator.

11) RRI_OscillLog_x_PVC_Freq (orig: oscill_log_x_pvc → new: RRI_OscillLog_x_PVC_Freq)
- Definition: log-oscillation × PVC frequency.
- Meaning: interaction of two strong single features.

12) RRI_Oscill_x_RecoveryCV (orig: oscill_x_cv → new: RRI_Oscill_x_RecoveryCV)
- Definition: oscillation index × CV of recovery half-life (feature 9).
- Meaning: synergy between oscillation and recovery instability.

13) PVC_Burden_NormIndex (orig: pvc_burden_index → new: PVC_Burden_NormIndex)
- Definition: normalized PVC burden: `(PVCs_per_hour/50) * (oscill_rr_median/0.5)` with cap at 2.0.
- Meaning: interpretable composite of frequency × oscillation.

14) RRI_Recovery_CapacityIndex (orig: recovery_capacity → new: RRI_Recovery_CapacityIndex)
- Definition: recovery capacity index: `fast_recovery_ratio / early_late_ratio`, where `fast_recovery_ratio = proportion(halflife≤10s)`.
- Meaning: larger means “more PVCs recover quickly + more balanced early/mid speeds”.

15) RRI_Oscill_Freq_LogRatio (orig: oscill_pvc_log_ratio → new: RRI_Oscill_Freq_LogRatio)
- Definition: log-oscillation / log PVC frequency: `oscill_log / log(max(PVCs_per_hour,1))`.
- Meaning: relative strength of oscillation vs frequency on the log scale.

Notes:

- All features are record-level aggregates; PVC-level metrics (half-life, oscillation, etc.) are computed per PVC, then aggregated by median/mean/ratios/interactions.
- Other computed but excluded “low-value features” (e.g., `halflife_rr_median`, `HRT_TS_low_frac`, `halflife_qtc_median`, `PVC_to_next_nonPVC_sec_*`) are fully listed in the script’s `lowValueFeatures`. The generator ultimately keeps only the 15 features above for training.
- Before training: group-wise stratified split by patient, missing imputation with training-set medians, MinMax normalization (fit on training), and optional demographics (disabled by default).

---

### Full Candidate Set and Exclusions (Transparency)

Before selecting the 15 core features, v3 constructs a comprehensive record-level candidate set from PVC recovery, HRT, QTc/T recovery, and complexity metrics, then applies transformations/ratios/interactions/CVs/threshold flags. Below are the main categories/examples and the explicit low-value exclusion list.

- Basic burden/quality:
	- PVCs_per_hour (per-hour PVC), PVC_count (total), recovery_failure_ratio (fraction not recovered within the max observation window), RR_Pre_RMSSD_mean (baseline RR RMSSD before PVC)
- Recovery time/rate:
	- halflife_rr_median / halflife_rr30_median / halflife_rr10_median (50%/30%/10% half-life medians)
	- tau_rr_median (time constant median), recovery_rate_hmean (harmonic mean of half-life as “speed”)
- Recovery shape/oscillation/stability:
	- oscill_rr_median (oscillation index median), late_var_rr (late variance)
	- PVC_to_next_nonPVC_sec_mean/max/min (time from PVC to the next non-PVC beat)
- HRT (heart rate turbulence):
	- HRT_TO_abnormal_frac (abnormal TO fraction), HRT_TS_low_frac (low TS fraction)
- T/QTc recovery:
	- halflife_tamp_median, halflife_qtc_median, auc_qtc_mean, qtc_overshoot_frac_mean, qtc_over_mag
- Complexity/nonlinearity/correlation:
	- sampen_rr (sample entropy), lzc_rr (Lempel-Ziv on binary first-difference signs), poincare_ratio_pp (SD1/SD2 ratio pre-post)
	- tamp_rr_corr_med (correlation between T amplitude and RR deviation)
- Slopes/segments:
	- rr_dev_slope_0_5s_med (0–5 s slope median), slope_rr_dev_5_15s (5–15 s slope)
- PVC coupling/HR dynamics:
	- coupling_ratio_mean (RR_pre / baseline RR), hr_peak_accel_bpm_mean (peak HR acceleration after PVC)
- Derived/interaction/transforms:
	- log/sqrt: oscill_log, pvc_freq_log, pvc_freq_sqrt, halflife30_log, halflife10_sqrt
	- interactions: oscill_x_pvc_freq, oscill_x_recovery, pvc_x_variability, oscill_log_x_pvc, oscill_x_cv, early_late_x_cv
	- ratios/normalized: oscill_per_pvc, recovery_efficiency, early_late_ratio, oscill_pvc_log_ratio
	- variability: oscill_cv, halflife_cv, combined_cv, recovery_oscill_cv_ratio, tau_x_cv
	- composite/threshold: composite_risk_score, instability_score, pvc_burden_index, recovery_capacity, high_freq_high_oscill, slow_recovery_high_cv, composite_high_risk, HRT_abnormal_combined, HRT_risk_category

To improve generalization and stability, the following low-value features are excluded before final training (from `lowValueFeatures` with SHAP/AUC/multivariate checks):

- Batch 1 (native low-value/redundant):
	- PVC_count, coupling_ratio_mean, recovery_rate_hmean, halflife_rr_median, PVC_to_next_nonPVC_sec_mean, PVC_to_next_nonPVC_sec_min, HRT_TS_low_frac, auc_qtc_mean, qtc_overshoot_frac_mean, halflife_tamp_median, halflife_qtc_median, rr_dev_late_minus_early
- Batch 2 (newly introduced but poor, e.g., SHAP<0.015):
	- halflife_rr10_median, tamp_rr_corr_med, halflife10_sqrt, recovery_efficiency, oscill_per_pvc, recovery_deceleration, HRT_risk_category, pvc_x_variability
- Batch 3 (not in Top15 or likely overfitting from complex interactions):
	- oscill_x_recovery, halflife30_log, combined_cv, early_late_x_cv, instability_score, high_freq_high_oscill, slow_recovery_high_cv, composite_high_risk, RR_Pre_RMSSD_mean, recovery_failure_ratio, halflife_rr30_median
- Batch 4 (final pruning to 15):
	- tau_x_cv, HRT_abnormal_combined, pvc_freq_log, rr_dev_slope_0_5s_med, oscill_cv, hr_peak_accel_bpm_mean, recovery_oscill_cv_ratio, PVC_to_next_nonPVC_sec_max, fast_recovery_ratio

The 15 retained training features are fully explained in the previous section (with renamed names).

---

## Training and Evaluation Setup (v3)

- Algorithm: AdaBoostM1 (auto-search indicates best among candidates; also considers `LogitBoost/GentleBoost/RUSBoost`)
- Typical hyperparameters: NLC=250, LR=0.015, MaxSplits=20, MinLeaf=60, NumVars='sqrt'
- Class imbalance: cost matrix (FP=2, FN=8) to emphasize recall for the Dead class
- Threshold moving: optimize F1 (Dead) on cross-validation; probabilities may use Platt calibration
- Cross-validation: K=10
- Metrics: AUC/AUPRC/Brier, confusion matrices (Dead as positive), calibration/PR/ROC/decision curves
- Importance: model gain, SHAP (test subset sampling, MATLAB R2021a+ or local implementation), univariate discriminative power (AUC/MI/Cohen’s d/rank-sum)

Main outputs:

- `sca_classifier_v3.mat` contains: model, required feature names, key train/test metrics, univariate effectiveness, importances
- `roc_curve_train_test.png`, `pr_curve_train_test.png`, `calibration_curve_train_test.png`, `decision_curve_test.png`
- `feature_effectiveness_*.csv`, `shap_importance_v3.csv`, `ablation_demographic_physio_combined_v3.csv`

---

## Suggested Repro Steps

In MATLAB:

1) Ensure data are placed as described under “Data layout requirement”, and set the project root (`pwd`) to this repo root.
2) Run `generate_beats_classifier.m` (if you need to retrain the beat-type model; requires rpoints).
3) Run `predict_shhs1_all.m` to produce `*_info.mat` (no rpoints needed).
4) Run `generate_sca_classifier_v3.m` to produce features, train, and evaluate.

Visualization: `view_shhs_ecg_v2.m` / `view_shhs_ecg_v3.m` help inspect records and debug.

---

## Important Implementation Details and Edge Cases

- Baseline/recovery statistics: within `baselineSec=50s`, take non-PVC and good-SQI beats to compute means/SDs for RR/T/QTc; fallback to global stats when insufficient.
- Recovery determination: dev sequence entering thresholds (e.g., 50%/30%/10%) for `consecStableBeats=10` consecutive beats is considered reaching the corresponding “half-life”.
- Reasonable RR range: `[0.30, 2.50] s`; PVC quality gating at record level: `minPVCPerRecord=10` and `minSQIRatio=0.60`.
- QTc approximation: use global T indices and QRS duration to estimate QT; Bazett cube-root correction; range capping; derive deviation and over-threshold fractions.
- Group-wise stratified split: by patient ID to avoid data leakage across train/test.

---

## Dependencies and Environment

- MATLAB R2021a+ (recommended)
- Toolboxes: Statistics and Machine Learning, Signal Processing
- (Optional) WFDB MATLAB toolbox in `mcode/` (bundled copy; mainly for extension/validation; not a hard dependency for v3)

---

## FAQ

- Q: No rpoints available—can I still run v3?
	- A: Yes. Use an existing `results/trainedClassifier_latest.mat` (or train your own) and run `predict_shhs1_all.m`, which needs no rpoints and will classify PVC and generate `*_info.mat`. Then run `generate_sca_classifier_v3.m`.

- Q: Feature names mismatch?
	- A: v3 renames the 15 core features before training (for better statistical/physiological semantics). The script contains the original→new mapping, and the list above shows both.

- Q: Why were some QTc/T-wave related features removed?
	- A: Based on SHAP/univariate power and multivariate validation, low-performing or redundant features were removed to keep 15 core features for better generalization.

---

## Acknowledgments

- SHHS data from the National Sleep Research Resource (NSRR)
- WFDB tools from PhysioNet (see `mcode/` and its LICENSE)

---

## Changelog (aligned with scripts)

- 2025-09-04: v3 initial
- 2025-10-08: refactor and comments; centralized config; consistent feature selection and naming

For more details or script parameters, please refer to headers and inline comments in the corresponding MATLAB files.

