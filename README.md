# ENG 7002 ECG Classification & SCA Prediction

This project provides a complete MATLAB pipeline for SHHS1 ECG signal processing, including filtering, Q/R/S/T delineation, beat quality (SQI) assessment, feature extraction, PVC vs. Other beat classification, and sudden cardiac arrest (SCA) risk–related post‑PVC window aggregation and modeling. Interactive visualization GUIs are included.

## 1. Top-Level Code & Directory Structure (Text View)

```
ENG-7002-ECG-Classification-SCA-Prediction/
├─ README.md                        ← Project documentation (this file)
├─ assessBeatsQuality.m             ← Beat template SQI evaluation (template build + correlation threshold gating)
├─ detectAndClassifyHeartbeats.m    ← Main Q/R/S/T detection + annotation matching + SQI filtering (calls assessBeatsQuality)
├─ ecgFilter.m                      ← ECG preprocessing (notch / baseline removal / low-pass / smoothing / QRS enhancement)
├─ extractHeartbeatFeatures.m       ← Derive RR / amplitude / duration / area features from beatInfo
├─ trainBeatsClassifier.m           ← Train PVC vs Other beat classifier (AdaBoost)
├─ generate_beats_classifier.m      ← Script: batch feature extraction + train/evaluate/save beat classifier
├─ predict_shhs1_all.m              ← Script: run detection + features + prediction for all EDFs using trained model
├─ generate_sca_classifier.m        ← Script: aggregate 5 s post-PVC features from *_info.mat (predPVCIndices, etc.) and train SCA model
├─ view_shhs_ecg_v1.m               ← GUI: raw / filtered ECG + annotated R peaks
├─ view_shhs_ecg_v2.m               ← GUI: filtered ECG + Q/R/S/T markers (internal detection); bottom raw + annotations
├─ view_shhs_ecg_v3.m               ← GUI: integrated viewer + on‑the‑fly feature extraction + model inference (PVC / risk)
├─ generate_sca_classifier.m        ← SCA classifier training script (extensive window feature engineering)
├─ mcode/                           ← External / auxiliary WFDB & HRV toolkit (PhysioNet MATLAB tools + examples)
│  ├─ *.m / *.jar / example/ ...    ← rdann, rdsamp, gqrs, ecgpuwave etc. (not directly invoked in core pipeline yet)
│  └─ html/                         ← Corresponding function HTML help
└─ results/                         ← Generated output (models & feature tables, created at runtime)
```

## 2. Main Data Flow & Call Relationships (Bottom → Top)

### 2.1 Low-Level Signal Processing & Detection
```
ecgFilter
	├─ (internal) apply_FIR_filter / apply_IIR_filter / apply_enhanced_pipeline
	└─ Produces filtered ECG (optionally QRS‑enhanced signal)

detectAndClassifyHeartbeats (core beat detection + structured beatInfo)
	├─ local_detect_r_hybrid_fast (hybrid R strategy)
	├─ local_remove_t_as_r (remove T misdetections)
	├─ local_detect_qs (Q/S + QRS boundaries)
	├─ local_detect_twave (T tagging)
	├─ assessBeatsQuality (template SQI gating)
	└─ Output: heartbeatSegments, beatInfo (segmentStartIndex / rIndex / qIndex / sIndex / tIndex / qrsOn/Off)

assessBeatsQuality (beat segment quality scoring)
	├─ Build templates: Other / PVC / Global (subject to minimum counts)
	├─ maxCorrWithShift sliding alignment correlation
	└─ Output: isGood mask + corrValues + template class used

extractHeartbeatFeatures
	├─ Input: beatInfo + fs
	├─ Compute: RR intervals, Q/R/S amplitudes, QRS duration, QRS area (prefer qrsOn/Off)
	└─ Output: featureTable (includes BeatType)
```

### 2.2 Beat Classification (PVC vs Other)
```
generate_beats_classifier (script)
	├─ Enumerate EDF + rpoint annotations
	├─ ecgFilter → detectAndClassifyHeartbeats → extractHeartbeatFeatures
	├─ Merge all features → stratified split (train/test)
	├─ trainBeatsClassifier (AdaBoost with optional cost & threshold moving)
	├─ Evaluate: confusion matrix / metrics
	└─ Save: trainingFeatureTable.mat / testingFeatureTable.mat / trainedClassifier_latest.mat

trainBeatsClassifier
	├─ Class frequency weights + optional cost matrix
	├─ fitcensemble (AdaBoostM1 + templateTree)
	├─ Optional threshold-adjusted predictFcn wrapper
	└─ Output: trainedClassifier struct + 10-fold CV accuracy
```

### 2.3 Batch Prediction & Persistence
```
predict_shhs1_all (script)
	├─ Load trainedClassifier_latest.mat
	├─ For each EDF: ecgFilter → detectAndClassifyHeartbeats (no annotations) → extractHeartbeatFeatures
	├─ Median-impute missing numeric features → model predict (optional external PVC threshold override)
	└─ Save per-record *_beats_features_pred.mat (allBeatInfo / raw & imputed features / predicted labels / stats)
```

### 2.4 SCA Risk Feature Aggregation & Classification
```
generate_sca_classifier (script)
	├─ Depends on *_info.mat (predPVCIndices, patientVital)
	├─ Read EDF → ecgFilter → detectAndClassifyHeartbeats for global R series
	├─ For each PVC build post-PVC 5 s window features (stat / rhythm / morphology)
	├─ Aggregate feature table → train (local_train_sca_classifier: AdaBoost + weights + cost + 5-fold CV)
	├─ Threshold moving / evaluation (accuracy, sensitivity, specificity, AUC)
	└─ Save post_ectopic_features.mat / sca_classifier_post_ectopic.mat
```

### 2.5 Visualization GUIs
```
view_shhs_ecg_v1
	├─ Open EDF → ecgFilter (method 2, skip notch) → plot raw & inverted filtered signals
	├─ Try load rpoint annotations (seconds / indices) → show R markers & labels

view_shhs_ecg_v2
	├─ Open EDF → ecgFilter → detectAndClassifyHeartbeats (no annotations) for Q/R/S/T
	├─ Top: detection signal with Q/R/S/T markers; bottom: raw + R annotations

view_shhs_ecg_v3
	├─ Advanced window / navigation / batch ops
	├─ Buttons trigger detectAndClassifyHeartbeats + extractHeartbeatFeatures
	├─ Load trainedClassifier_latest for online PVC classification
	└─ Can be extended with SCA risk model
```

## 3. End-to-End Typical Workflow
```
1. (Optional) Run generate_beats_classifier.m to train/update PVC model
2. Run predict_shhs1_all.m to produce *_beats_features_pred.mat for all records
3. Generate *_info.mat (external/preceding script producing predPVCIndices & patientVital)
4. Run generate_sca_classifier.m to train SCA risk model
5. Use view_shhs_ecg_v2 / v3 for interactive inspection & validation
```

## 4. Key Data Structures
```
beatInfo(i):
	.beatType          Beat label (PVC / Other; default Other if no annotations)
	.segment           Beat waveform segment
	.segmentStartIndex Segment start (global sample index)
	.rIndex/qIndex/sIndex/tIndex  Indices within segment
	.qrsOnIndex/.qrsOffIndex      QRS onset/offset (segment)
	.sqiIsGood / .sqiCorr / .sqiTemplateClass  SQI quality metrics

featureTable row:
	RR_Prev / RR_Post / R_Amplitude / Q_Amplitude / S_Amplitude / QRS_Duration / QRS_Area / BeatType
```

### 4.1 Beat-Level Features Used in `generate_beats_classifier.m`
Source: output of `extractHeartbeatFeatures` + BeatType label; `trainBeatsClassifier` defines `predictorNames`.
```
Predictor (numeric) columns:
	RR_Prev          Previous RR interval (s)
	RR_Post          Following RR interval (s)
	R_Amplitude      R peak amplitude (raw segment value)
	Q_Amplitude      Q amplitude
	S_Amplitude      S amplitude
	QRS_Duration     QRS duration (prefer qrsOn/off else fallback)
	QRS_Area         Absolute QRS area (simple integral)
Label:
	BeatType (PVC / Other)
```

### 4.2 Record-Level (Post-PVC 5 s) Feature Set in `generate_sca_classifier.m`
The script builds per-PVC window primitive metrics, then applies statistical aggregation (mean / std / median / min / max) and adds derived ratio / proportion / variability indicators. Final predictors are numeric columns excluding `record`, `pvc_r_index` (if present), and the response `patientVital` (selected by `local_training_columns`).

Primitive per-window features (each PVC window):
```
RR_Pre, RR_Post1, RR_Post2
HRT_TO = (RR_Post1 - RR_Pre)/RR_Pre        (Heart Rate Turbulence Onset approximation)
HRT_TS = (RR_Post2 - RR_Pre)/2             (Turbulence Slope approximation)
QRS_Dur_PVC, R_Amp_PVC
Beats_in_5s (R count within window)
HR_Pre = 60 / RR_Pre
HR_Post1 = 60 / RR_Post1
HR_5s = Beats_in_5s / windowSeconds * 60
HR_Accel = HR_Post1 - HR_Pre
CompRatio = RR_Post1 / RR_Pre
PVC_Interval (adjacent PVC interval series)
```

Statistical aggregation applied to each numeric series (except certain proportion / derived only fields):
```
<Name>_mean, <Name>_std, <Name>_median, <Name>_min, <Name>_max
Examples: RR_Pre_mean, RR_Pre_std, ..., HRT_TO_mean, QRS_Dur_PVC_max, HR_5s_median
```

Additional derived / proportion / variability metrics:
```
PVCs_per_hour                 PVC count per hour
HRT_TO_neg_frac               Fraction with HRT_TO < 0
QRS_Prolonged_frac            Fraction QRS_Dur_PVC > 0.12 s
RR_Pre_CV, RR_Post1_CV        RR coefficient of variation (std/mean)
RR_Pre_RMSSD, RR_Post1_RMSSD  RR RMSSD (short-term variability)
HR_Accel_*                    Heart rate acceleration stats (mean/std/median/min/max)
CompRatio_*                   Compensatory ratio stats
PVC_Interval_*                PVC interval stats
record_pvc_count              Total PVC count in record
patientVital                  Label (1=Alive, 0=Dead) — response variable
```
