% File: generate_sca_classifier.m
% Type: Script
% Description:
%   Using pre-generated *_info.mat (predPVCIndices and patientVital), aggregate post-PVC 5s features
%   per record, stratify, train and evaluate SCA risk classifier, and save model and features.
% Usage:
%   Run after generating *_info.mat with batch script.
% Inputs/Dependencies:
%   - Data dir: shhs/polysomnography/edfs/shhs1
%   - Functions: ecgFilter, detectAndClassifyHeartbeats
% Outputs:
%   - results/post_ectopic_features.mat
%   - results/sca_classifier_post_ectopic.mat
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

clc; clear; close all;

% ===================== Configuration =====================
rootDir = pwd;
edfDir = fullfile(rootDir, 'shhs','polysomnography','edfs','shhs1');
if ~isfolder(edfDir)
    error('EDF directory not found: %s', edfDir);
end
addpath(genpath(rootDir));

% SHHS1 sampling rate
fs = 125; 

% Process only first N EDFs with *_info.mat (<=0 means all)
processFirstN = 30;

% Post-PVC window length (seconds)
postWindowSec = 5.0;

% Preprocessing: filter entire record once (using ecgFilter)

% Training configuration
splitRatioTest = 0.20; % 80/20
randomSeed = 42;
enableCostMatrix = true; costFP = 1; costFN = 3; % heavier penalty for false negatives
enableThresholdMoving = true; pvcThreshold = 0.40; % probability threshold (Alive=1, Dead=0 → used for SCA risk)

% Outputs
resultsDir = fullfile(rootDir, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end
featOutFile = fullfile(resultsDir, 'post_ectopic_features.mat');
modelOutFile = fullfile(resultsDir, 'sca_classifier_post_ectopic.mat');

fprintf('=== Post-ectopic 5s feature extraction and SCA classifier training ===\n');
fprintf('EDF dir: %s\n', edfDir);

% ===================== Enumerate EDF and _info.mat =====================
edfFiles = dir(fullfile(edfDir, 'shhs1-*.edf'));
if isempty(edfFiles)
    error('No EDF files found: %s', edfDir);
end

% Keep only records with *_info.mat
edfPaths = {};
infoPaths = {};
for i = 1:numel(edfFiles)
    [~, nm, ~] = fileparts(edfFiles(i).name);
    infoPath = fullfile(edfDir, sprintf('%s_info.mat', nm));
    if exist(infoPath, 'file')
        edfPaths{end+1,1} = fullfile(edfDir, edfFiles(i).name); %#ok<AGROW>
        infoPaths{end+1,1} = infoPath; %#ok<AGROW>
    end
end

if isempty(edfPaths)
    error('No EDF with *_info.mat found. Please run predict_shhs1_rpoints.m first');
end

% Sort and optionally truncate
[edfPaths, sortIdx] = sort(edfPaths);
infoPaths = infoPaths(sortIdx);
if processFirstN > 0
    edfPaths = edfPaths(1:min(processFirstN, numel(edfPaths)));
    infoPaths = infoPaths(1:min(processFirstN, numel(infoPaths)));
end

fprintf('Will process %d records (with *_info.mat).\n', numel(edfPaths));

% ===================== Iterate records and generate samples =====================
allFeatureTables = cell(0,1);
numTotalPVCWindows = 0;

for iFile = 1:numel(edfPaths)
    edfPath = edfPaths{iFile};
    infoPath = infoPaths{iFile};
    [edfFolder, edfName, ~] = fileparts(edfPath);
    recBase = edfName; % shhs1-XXXXXX
    fprintf('\n=== [%d/%d] Processing %s ===\n', iFile, numel(edfPaths), [edfName '.edf']);

    % Load predPVCIndices and patientVital
    S = load(infoPath);
    if ~isfield(S,'predPVCIndices') || isempty(S.predPVCIndices)
        fprintf('  Skip: no predPVCIndices.\n');
        continue;
    end
    if ~isfield(S,'patientVital') || ~isfinite(S.patientVital)
        fprintf('  Skip: no patientVital.\n');
        continue;
    end
    predPVCIndices = double(S.predPVCIndices(:));
    patientVital = double(S.patientVital); % 0/1

    % Read EDF and extract ECG channel
    try
        TT = edfread(edfPath);
    catch ME
        fprintf('  Skip: edfread failed: %s\n', ME.message);
        continue;
    end
    varNames = TT.Properties.VariableNames;
    ecgIdx = find(strcmp(varNames, 'ECG'), 1);
    if isempty(ecgIdx)
        vlow = lower(varNames);
        ecgIdx = find(contains(vlow, 'ecg') | contains(vlow, 'ekg'), 1, 'first');
    end
    if isempty(ecgIdx)
        fprintf('  Skip: ECG channel not found.\n');
        continue;
    end
    ecgVarName = varNames{ecgIdx};
    ecgCol = TT.(ecgVarName);
    if iscell(ecgCol)
        try
            ecg = vertcat(ecgCol{:});
        catch
            ecg = [];
            for k = 1:numel(ecgCol)
                v = ecgCol{k};
                if isstring(v) || ischar(v)
                    v = str2num(v); %#ok<ST2NM>
                end
                ecg = [ecg; v(:)]; %#ok<AGROW>
            end
        end
    elseif isnumeric(ecgCol)
        ecg = ecgCol(:);
    else
        fprintf('  Skip: ECG channel type unsupported: %s\n', class(ecgCol));
        continue;
    end
    ecg = double(ecg);
    N = numel(ecg);
    fprintf('  ECG samples=%d, fs=%d Hz, PVC count=%d\n', N, fs, numel(predPVCIndices));

    % Filter entire record (using ecgFilter)
    try
        [ecgFiltered, ~] = ecgFilter(ecg, fs, 2, 0);
    catch ME
        fprintf('  Filtering warning: %s, using raw signal.\n', ME.message);
        ecgFiltered = ecg;
    end

    % Detect R/QRS across entire record (skip P/T detection)
    try
        ATRTIMED = [];
        ANNOTD = {};
        [~, beatInfo, stats] = detectAndClassifyHeartbeats(ecgFiltered, ATRTIMED, ANNOTD, fs); %#ok<ASGLU>
    catch ME
        fprintf('  Skip: beat detection error: %s\n', ME.message);
        continue;
    end
    if isempty(beatInfo)
        fprintf('  Skip: no valid beats.\n');
        continue;
    end
    % Global R indices (samples)
    rGlobalAll = zeros(numel(beatInfo),1);
    for kk2 = 1:numel(beatInfo)
        rGlobalAll(kk2) = double(beatInfo(kk2).segmentStartIndex) + double(beatInfo(kk2).rIndex) - 1;
    end

    % For each PVC, slice post-PVC 5-second windows (vectorized)
    winLen = round(postWindowSec * fs);
    recordPVCCount = numel(predPVCIndices);

    % Sort PVC indices for sliding window counting
    pvcIndicesSorted = sort(predPVCIndices(:));
    numPVC = numel(pvcIndicesSorted);
    if numPVC == 0
        fprintf('  No valid PVC window.\n');
        continue;
    end

    % Map PVC R indices to beatInfo indices (vectorized nearest neighbor)
    idxPVC_all_sorted = round(interp1(rGlobalAll, 1:numel(rGlobalAll), pvcIndicesSorted, 'nearest', 'extrap'));
    idxPVC_all_sorted = max(1, min(numel(rGlobalAll), idxPVC_all_sorted));

    % Precompute RR vectors
    numBeatsRec = numel(rGlobalAll);
    rr_pre_vec = nan(numBeatsRec,1); rr_post1_vec = nan(numBeatsRec,1); rr_post2_vec = nan(numBeatsRec,1);
    if numBeatsRec >= 2
        rr_pre_vec(2:end) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
        rr_post1_vec(1:end-1) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
    end
    if numBeatsRec >= 3
        rr_post2_vec(1:end-2) = (rGlobalAll(3:end) - rGlobalAll(1:end-2)) / fs;
    end

    % Precompute QRS duration and R amplitude
    qrs_dur_vec = nan(numBeatsRec,1);
    r_amp_vec = nan(numBeatsRec,1);
    for ii = 1:numBeatsRec
        b = beatInfo(ii);
        on = NaN; off = NaN;
        if isfield(b,'qrsOnIndex'), on = double(b.qrsOnIndex); end
        if isfield(b,'qrsOffIndex'), off = double(b.qrsOffIndex); end
        if isfinite(on) && isfinite(off) && off > on
            qrs_dur_vec(ii) = (off - on) / fs;
        end
        if isfield(b,'segment') && isfield(b,'rIndex') && ~isempty(b.segment) && isfinite(b.rIndex) && b.rIndex>=1 && b.rIndex<=numel(b.segment)
            r_amp_vec(ii) = double(b.segment(b.rIndex));
        end
    end

    % Calculate end point and duration for each PVC window
    winStartVec = pvcIndicesSorted;
    winEndVec = min(N, winStartVec + winLen);
    winDurSecVec = (winEndVec - winStartVec) / fs;

    % Use two pointers to count R waves within window (Beats_in_5s)
    counts = zeros(numPVC,1);
    left = 1; right = 0;
    for kk = 1:numPVC
        a = winStartVec(kk); b = winEndVec(kk);
        while right < numBeatsRec && rGlobalAll(right+1) <= b
            right = right + 1;
        end
        while left <= numBeatsRec && rGlobalAll(left) < a
            left = left + 1;
        end
        if right >= left
            counts(kk) = right - left + 1;
        else
            counts(kk) = 0;
        end
    end

    % Assemble features (vectorized indexing)
    RR_Pre = rr_pre_vec(idxPVC_all_sorted);
    RR_Post1 = rr_post1_vec(idxPVC_all_sorted);
    RR_Post2 = rr_post2_vec(idxPVC_all_sorted);
    HRT_TO = (RR_Post1 - RR_Pre) ./ (RR_Pre + eps);
    HRT_TS = (RR_Post2 - RR_Pre) / 2;
    QRS_Dur_PVC = qrs_dur_vec(idxPVC_all_sorted);
    R_Amp_PVC = r_amp_vec(idxPVC_all_sorted);
    Beats_in_5s = counts;

    % Heart rate features
    HR_Pre = 60 ./ RR_Pre;
    HR_Post1 = 60 ./ RR_Post1;
    HR_5s = (counts ./ max(winDurSecVec, eps)) * 60;

    % Record-level columns
    recordCol = repmat({recBase}, numPVC, 1);
    vitalCol = repmat(patientVital, numPVC, 1);
    pvcIndexCol = pvcIndicesSorted;
    recPVCcountCol = repmat(recordPVCCount, numPVC, 1);

    % ===================== Record-level aggregated features =====================
    % Additional derived quantities
    HR_Accel = HR_Post1 - HR_Pre;                  % Heart rate acceleration after PVC
    CompRatio = RR_Post1 ./ (RR_Pre + eps);        % Compensatory pause ratio
    pvcIntervals = diff(pvcIndicesSorted) / fs;    % PVC-PVC intervals (seconds)

    % SQI features removed (not included in record-level features as per requirements)

    % Record duration and PVC density
    recordDurationHr = (N / fs) / 3600;
    pvcPerHour = numPVC / max(recordDurationHr, eps);

    % Aggregate to struct (for easy extension)
    Srec = struct();
    Srec.record = recBase;
    Srec.patientVital = patientVital;
    Srec.record_pvc_count = recordPVCCount;
    Srec.PVCs_per_hour = pvcPerHour;

    % Perform statistical aggregation on core window-level metrics (mean/std/median/min/max)
    Srec = local_add_stats(Srec, RR_Pre,    'RR_Pre');
    Srec = local_add_stats(Srec, RR_Post1,  'RR_Post1');
    Srec = local_add_stats(Srec, RR_Post2,  'RR_Post2');
    Srec = local_add_stats(Srec, HRT_TO,    'HRT_TO');
    Srec = local_add_stats(Srec, HRT_TS,    'HRT_TS');
    Srec = local_add_stats(Srec, QRS_Dur_PVC, 'QRS_Dur_PVC');
    Srec = local_add_stats(Srec, R_Amp_PVC, 'R_Amp_PVC');
    Srec = local_add_stats(Srec, Beats_in_5s, 'Beats_in_5s');
    Srec = local_add_stats(Srec, HR_Pre,    'HR_Pre');
    Srec = local_add_stats(Srec, HR_Post1,  'HR_Post1');
    Srec = local_add_stats(Srec, HR_5s,     'HR_5s');

    % SCA-related additional aggregated features
    % - HRT-TO negative value proportion (healthy usually more negative; non-negative ↑ may indicate ↑ risk)
    mask_to = isfinite(HRT_TO);
    if any(mask_to)
        Srec.HRT_TO_neg_frac = mean(double(HRT_TO(mask_to) < 0));
    else
        Srec.HRT_TO_neg_frac = NaN;
    end
    % - QRS prolongation proportion (>120 ms)
    mask_qrs = isfinite(QRS_Dur_PVC);
    if any(mask_qrs)
        Srec.QRS_Prolonged_frac = mean(double(QRS_Dur_PVC(mask_qrs) > 0.12));
    else
        Srec.QRS_Prolonged_frac = NaN;
    end
    % - Heart rate acceleration/compensation ratio
    Srec = local_add_stats(Srec, HR_Accel,   'HR_Accel');
    Srec = local_add_stats(Srec, CompRatio,  'CompRatio');
    % - RR coefficient of variation (rhythm instability)
    Srec.RR_Pre_CV   = local_cv(RR_Pre);
    Srec.RR_Post1_CV = local_cv(RR_Post1);
    % - RR short-term variability RMSSD
    Srec.RR_Pre_RMSSD   = local_rmssd(RR_Pre);
    Srec.RR_Post1_RMSSD = local_rmssd(RR_Post1);
    % - PVC-PVC interval statistics
    Srec = local_add_stats(Srec, pvcIntervals, 'PVC_Interval');
    % - Signal quality (template-related SQI) removed as per requirements

    % Convert to single row table and accumulate
    TrecAgg = struct2table(Srec);
    allFeatureTables{end+1,1} = TrecAgg; %#ok<AGROW>
    numTotalPVCWindows = numTotalPVCWindows + numPVC;
    fprintf('  Valid PVC windows: %d\n', numPVC);
end

if isempty(allFeatureTables)
    error('No post-ectopic window features obtained.');
end

% Merge all
featureTable = vertcat(allFeatureTables{:});
fprintf('\nTotal samples: %d (each row per record)\n', height(featureTable));

% Save features
try
    save(featOutFile, 'featureTable', 'fs', 'postWindowSec');
    fprintf('✓ Saved features: %s\n', featOutFile);
catch ME
    fprintf('Failed to save features: %s\n', ME.message);
end

% ===================== Training/Evaluation =====================

% Construct training table: select numerical feature columns and response column
[predictorNames, responseName] = local_training_columns(featureTable);

% Simple missing value imputation (by column median)
for i = 1:numel(predictorNames)
    col = featureTable.(predictorNames{i});
    if ~isnumeric(col), continue; end
    nanMask = isnan(col);
    if any(nanMask)
        medv = median(col(~nanMask));
        if isempty(medv) || isnan(medv), medv = 0; end
        col(nanMask) = medv;
        featureTable.(predictorNames{i}) = col;
    end
end

% Stratified split (by patientVital)
rng(randomSeed);
y = featureTable.(responseName);
if iscell(y), y = cell2mat(y); end
if islogical(y), y = double(y); end
idxAlive = find(y == 1);
idxDead  = find(y == 0);
numAliveTest = round(numel(idxAlive) * splitRatioTest);
numDeadTest  = round(numel(idxDead)  * splitRatioTest);
idxAlivePerm = idxAlive(randperm(numel(idxAlive)));
idxDeadPerm  = idxDead(randperm(numel(idxDead)));
testIdx = [idxAlivePerm(1:min(numAliveTest,numel(idxAlivePerm))); idxDeadPerm(1:min(numDeadTest,numel(idxDeadPerm)))];
testMask = false(height(featureTable),1); testMask(testIdx) = true;
trainMask = ~testMask;

trainTbl = featureTable(trainMask, :);
testTbl  = featureTable(testMask, :);

fprintf('Training samples=%d, Test samples=%d\n', height(trainTbl), height(testTbl));

% Construct model options and train (reuse trainClassifier style: weighted + cost + threshold)
modelOptions = struct();
modelOptions.enableCostMatrix = logical(enableCostMatrix);
modelOptions.costFP = costFP;
modelOptions.costFN = costFN;
modelOptions.enableThresholdMoving = logical(enableThresholdMoving);
modelOptions.pvcThreshold = pvcThreshold;

[trainedClassifier, valAcc] = local_train_sca_classifier(trainTbl, predictorNames, responseName, modelOptions);
fprintf('Cross-validation accuracy (internal training): %.2f%%\n', 100*valAcc);

% Test set evaluation
[predLabel, scores] = local_predict_with_threshold(trainedClassifier, testTbl, predictorNames, responseName); %#ok<ASGLU>
yTrue = testTbl.(responseName);
if iscell(yTrue), yTrue = cell2mat(yTrue); end

% Build confusion matrix with Death(0) as positive class
yTrueDead = 1 - yTrue;      % 1=Dead,0=Alive
yHatDead  = 1 - predLabel;  % 1=Pred Dead,0=Pred Alive
TP = sum((yTrueDead==1) & (yHatDead==1));
FN = sum((yTrueDead==1) & (yHatDead==0));
FP = sum((yTrueDead==0) & (yHatDead==1));
TN = sum((yTrueDead==0) & (yHatDead==0));

% Text confusion matrix (rows: Actual, cols: Pred; order: Alive(1), Dead(0))
fprintf('\n— Text Confusion Matrix (positive=Death) —\n');
hdr = sprintf('%-12s%12s%12s', 'Actual\\Pred', 'Alive(1)', 'Dead(0)');
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, length(hdr)));
% Actual Alive(1): Pred Alive=TN, Pred Dead=FP
rowAlive = sprintf('%-12s%12d%12d', 'Alive(1)', TN, FP);
fprintf('%s\n', rowAlive);
% Actual Dead(0): Pred Alive=FN, Pred Dead=TP
rowDead  = sprintf('%-12s%12d%12d', 'Dead(0)', FN, TP);
fprintf('%s\n', rowDead);

% Metrics (positive=Death)
[acc, sens, spec, auc] = local_binary_metrics(yTrue, predLabel, scores);
prec = TP / max(TP + FP, eps);
if (prec + sens) > 0
    f1 = 2 * prec * sens / (prec + sens);
else
    f1 = NaN;
end
fprintf('\n— Performance (positive=Death) —\n');
fprintf('Accuracy: %.2f%%\n', 100*acc);
fprintf('Sensitivity (Death): %.2f%%\n', 100*sens);
fprintf('Specificity (Death): %.2f%%\n', 100*spec);
fprintf('Precision (Death): %.2f%%\n', 100*prec);
fprintf('F1 (Death): %.3f\n', f1);
fprintf('AUC(ROC): %.3f\n', auc);

% Save model
try
    save(modelOutFile, 'trainedClassifier', 'predictorNames', 'responseName', 'modelOptions');
    fprintf('✓ Saved model: %s\n', modelOutFile);
catch ME
    fprintf('Failed to save model: %s\n', ME.message);
end

fprintf('\nAll done.\n');

% ===================== Internal functions in this file =====================

function feat = local_extract_postPVC_features(rIdxGlobal, winStart, winEnd, fs, rGlobalAll, beatInfo)
% Extract lightweight post-ectopic 5-second window features (high relevance + high efficiency):
% - RR_pre: RR interval before PVC
% - RR_post1: RR interval from PVC to next R
% - RR_post2: RR interval from PVC to second R (if exists)
% - HRT-TO: (RR_post1 - RR_pre)/RR_pre
% - HRT-TS(approximate): (RR_post2 - RR_pre)/2
% - QRS_Dur_PVC, R_Amp_PVC (from corresponding beat in beatInfo)
% - Beats_in_5s: Number of R waves within window

    feat = struct();

    % Find PVC position in global R sequence
    [~, idxPVC] = min(abs(rGlobalAll - rIdxGlobal));
    % Previous and next R waves
    prevR = NaN; nextR1 = NaN; nextR2 = NaN;
    if ~isempty(idxPVC) && ~isnan(idxPVC)
        if idxPVC-1 >= 1, prevR = rGlobalAll(idxPVC-1); end
        if idxPVC+1 <= numel(rGlobalAll), nextR1 = rGlobalAll(idxPVC+1); end
        if idxPVC+2 <= numel(rGlobalAll), nextR2 = rGlobalAll(idxPVC+2); end
    end

    % RR features
    if isfinite(prevR)
        rr_pre = (rIdxGlobal - prevR)/fs;
    else
        rr_pre = NaN;
    end
    if isfinite(nextR1)
        rr_post1 = (nextR1 - rIdxGlobal)/fs;
    else
        rr_post1 = NaN;
    end
    if isfinite(nextR2)
        rr_post2 = (nextR2 - rIdxGlobal)/fs;
    else
        rr_post2 = NaN;
    end
    feat.RR_Pre = rr_pre;
    feat.RR_Post1 = rr_post1;
    feat.RR_Post2 = rr_post2;

    % HRT approximation
    if isfinite(rr_pre) && isfinite(rr_post1)
        feat.HRT_TO = (rr_post1 - rr_pre) / (rr_pre + eps);
    else
        feat.HRT_TO = NaN;
    end
    if isfinite(rr_pre) && isfinite(rr_post2)
        feat.HRT_TS = (rr_post2 - rr_pre)/2; % approximation
    else
        feat.HRT_TS = NaN;
    end

    % QRS duration and R amplitude (from corresponding beat in beatInfo, matched by closest R)
    % Index of corresponding beat in beatInfo
    [~, idxBeat] = min(abs(rIdxGlobal - (arrayfun(@(b) double(b.segmentStartIndex) + double(b.rIndex) - 1, beatInfo))));
    qrsDur = NaN; rAmp = NaN;
    if ~isempty(idxBeat) && ~isnan(idxBeat)
        b = beatInfo(idxBeat);
        if isfield(b,'qrsOnIndex') && isfield(b,'qrsOffIndex') && isfinite(b.qrsOnIndex) && isfinite(b.qrsOffIndex)
            if b.qrsOffIndex > b.qrsOnIndex
                qrsDur = (double(b.qrsOffIndex) - double(b.qrsOnIndex))/fs;
            end
        end
        if isfield(b,'segment') && isfield(b,'rIndex') && ~isempty(b.segment) && isfinite(b.rIndex) && b.rIndex>=1 && b.rIndex<=numel(b.segment)
            rAmp = double(b.segment(b.rIndex));
        end
    end
    feat.QRS_Dur_PVC = qrsDur;
    feat.R_Amp_PVC = rAmp;

    % Number of R waves within window (using global R sequence)
    feat.Beats_in_5s = sum(rGlobalAll >= winStart & rGlobalAll <= winEnd);
end

% Old slow auxiliary functions removed

function [predictorNames, responseName] = local_training_columns(T)
% Select columns for training:
% - predictor: numerical feature columns extracted above (exclude record/pvc_r_index/label fields)
% - response: patientVital (0/1)
    names = T.Properties.VariableNames;
    excl = ismember(names, {'record','pvc_r_index'});
    isNum = false(size(names));
    for i = 1:numel(names)
        isNum(i) = isnumeric(T.(names{i})) && isvector(T.(names{i}));
    end
    predictorMask = ~excl & isNum & ~strcmp(names, 'patientVital');
    predictorNames = names(predictorMask);
    responseName = 'patientVital';
end

function [trainedClassifier, validationAccuracy] = local_train_sca_classifier(trainingTable, predictorNames, responseName, modelOptions)
% Cost-sensitive + class weights + AdaBoost decision tree (consistent with project style)
    predictors = trainingTable(:, predictorNames);
    response = trainingTable.(responseName);
    if iscell(response), response = cell2mat(response); end
    % Observation weights (by class frequency)
    numObs = height(trainingTable);
    cntAlive = sum(response == 1);
    cntDead  = sum(response == 0);
    if cntAlive > 0 && cntDead > 0
        classCounts = [cntDead, cntAlive]; % order: 0,1
        classWeights = numObs ./ (2 .* classCounts);
        obsWeights = zeros(numObs,1);
        obsWeights(response == 0) = classWeights(1);
        obsWeights(response == 1) = classWeights(2);
    else
        obsWeights = ones(numObs,1);
    end

    template = templateTree('MaxNumSplits', 20, 'NumVariablesToSample', 'all');
    fitArgs = {predictors, response, 'Method','AdaBoostM1','NumLearningCycles',30,'Learners',template,'LearnRate',0.1,'Weights',obsWeights};
    if isfield(modelOptions,'enableCostMatrix') && modelOptions.enableCostMatrix
        % Cost matrix order by true class [0,1]
        C = [0, modelOptions.costFP; modelOptions.costFN, 0];
        fitArgs = [fitArgs, {'Cost', C}]; %#ok<AGROW>
    end
    cls = fitcensemble(fitArgs{:});

    % Simple cross-validation accuracy
    try
        cv = crossval(cls, 'KFold', 5);
        validationAccuracy = 1 - kfoldLoss(cv);
    catch
        validationAccuracy = NaN;
    end

    trainedClassifier = struct();
    trainedClassifier.ClassificationEnsemble = cls;
    trainedClassifier.RequiredVariables = predictorNames;
    trainedClassifier.Options = modelOptions;
end

function [predLabel, scores] = local_predict_with_threshold(model, T, predictorNames, ~)
    X = T(:, predictorNames);
    try
        [~, rawScores] = predict(model.ClassificationEnsemble, X);
    catch
        [~, rawScores] = predict(model.ClassificationEnsemble, X{:,:});
    end
    % fitcensemble returns two columns of scores for binary classification, column 2 usually corresponds to positive class (higher label value).
    % We treat patientVital==1 as positive class (alive), but risk prediction focuses more on patientVital==0.
    % For threshold moving convenience, define scores = P(y==0) (death risk).
    if size(rawScores,2) >= 2
        % MATLAB arranges columns by class order, need to extract class labels
        try
            cls = model.ClassificationEnsemble.ClassNames;
        catch
            cls = [0;1];
        end
        % Find column representing class 0
        col0 = find(cls==0, 1, 'first');
        if isempty(col0), col0 = 1; end
        scores = rawScores(:, col0);
    else
        scores = rawScores(:,1);
    end
    thr = 0.50;
    if isfield(model, 'Options') && isfield(model.Options, 'enableThresholdMoving') && model.Options.enableThresholdMoving
        if isfield(model.Options,'pvcThreshold') && isfinite(model.Options.pvcThreshold)
            thr = model.Options.pvcThreshold;
        end
    end
    predLabel = double(scores >= thr); % 1 indicates predicted risk (death); will flip later to align with yTrue (0/1=Alive/Dead)
    % Flip to patientVital semantics: we want to output yhat_vital (1=Alive, 0=Dead)
    predLabel = 1 - predLabel;
end

function [acc, sens, spec, auc] = local_binary_metrics(yTrueVital, yHatVital, scoresForDead)
% y*: 1=Alive, 0=Dead
    yTrueVital = double(yTrueVital(:));
    yHatVital = double(yHatVital(:));
    acc = mean(yTrueVital == yHatVital);
    % Sensitivity/specificity for "Dead" class
    yTrueDead = 1 - yTrueVital;
    yHatDead  = 1 - yHatVital;
    TP = sum((yTrueDead==1) & (yHatDead==1));
    FN = sum((yTrueDead==1) & (yHatDead==0));
    FP = sum((yTrueDead==0) & (yHatDead==1));
    TN = sum((yTrueDead==0) & (yHatDead==0));
    if TP+FN>0, sens = TP/(TP+FN); else, sens = NaN; end
    if TN+FP>0, spec = TN/(TN+FP); else, spec = NaN; end
    % AUC (based on scoresForDead)
    try
        [~,~,~,auc] = perfcurve(yTrueDead, scoresForDead, 1);
    catch
        auc = NaN;
    end
end


function S = local_add_stats(S, v, baseName)
% Append statistics of vector v to struct S, with field names like baseName_mean/std/median/min/max
    if nargin < 3 || isempty(baseName)
        baseName = 'feat';
    end
    if isempty(v)
        vals = [NaN, NaN, NaN, NaN, NaN];
    else
        v = v(:);
        v = v(isfinite(v));
        if isempty(v)
            vals = [NaN, NaN, NaN, NaN, NaN];
        else
            m = mean(v);
            s = std(v);
            med = median(v);
            mn = min(v);
            mx = max(v);
            vals = [m, s, med, mn, mx];
        end
    end
    S.([baseName '_mean']) = vals(1);
    S.([baseName '_std']) = vals(2);
    S.([baseName '_median']) = vals(3);
    S.([baseName '_min']) = vals(4);
    S.([baseName '_max']) = vals(5);
end

function cv = local_cv(v)
% Coefficient of variation: std/|mean| (ignore NaN). Return NaN if insufficient samples or mean too small.
    if isempty(v)
        cv = NaN; return; end
    v = v(:);
    v = v(isfinite(v));
    if numel(v) < 2
        cv = NaN; return; end
    mu = mean(v);
    sigma = std(v);
    den = abs(mu);
    if den < 1e-12
        cv = NaN; return; end
    cv = sigma / den;
end

function rmssd = local_rmssd(v)
% RR short-term variability: square root of mean of squared successive differences.
    if isempty(v) || numel(v) < 2
        rmssd = NaN; return; end
    v = v(:);
    v = v(isfinite(v));
    if numel(v) < 2
        rmssd = NaN; return; end
    d = diff(v);
    rmssd = sqrt(mean(d.^2));
end


