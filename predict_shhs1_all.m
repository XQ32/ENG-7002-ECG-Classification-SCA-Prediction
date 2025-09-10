% ========================================================================
% File: predict_shhs1_all.m
% Overview: Iterate over all SHHS1 EDFs (no rpoints labels) to detect beats
%           and predict beat types, and export each record's predicted PVC
%           global R indices and patient vital into *_info.mat.
% Responsibilities:
%   1) Enumerate all .edf files under shhs/polysomnography/edfs/shhs1.
%   2) Load the trained beat-type classifier (results/trainedClassifier_latest.mat).
%   3) Read the single-lead ECG channel (exact 'ECG' → fuzzy-match fallback)
%      and flatten to a double column vector.
%   4) Preprocess: ecgFilter(method=2, power_line_freq=0) (SHHS already has 60 Hz notch).
%   5) Beat detection detectAndClassifyHeartbeats (no external labels) to obtain beatInfo.
%   6) Feature extraction extractHeartbeatFeatures → impute missing values by column medians.
%   7) Use the model (optional threshold override pvcThresholdOverride) to predict each beat type ('PVC'/'Other').
%   8) Aggregate and save predPVCIndices (global R sample indices predicted as PVC) and patientVital into *_info.mat in the same directory.
% Inputs/Dependencies:
%   - Directory: shhs/datasets/shhs-cvd-summary-dataset-0.21.0.csv (for mapping vital: 0=Dead, 1=Alive)
%   - results/trainedClassifier_latest.mat (contains trainedClassifier or trainedModelPackage)
%   - Functions: ecgFilter, detectAndClassifyHeartbeats, extractHeartbeatFeatures
% Outputs:
%   - For each EDF: {record}_info.mat (predPVCIndices, patientVital)
% Tunable internal parameters:
%   fs (fixed 125 Hz), pvcThresholdOverride (override model PVC threshold; higher → higher precision, lower recall)
% Key implementation notes:
%   - Skip records if required model variables are missing (robustness).
%   - Compatible with two model save structures (trainedModelPackage / direct trainedClassifier).
%   - Vital info mapped via nsrrid→vital (patientVital=NaN if missing).
%   - Before prediction, median-impute required feature columns to avoid NaNs breaking prediction.
% Edge cases and robustness:
%   - If any of edfread/filter/detect/feature/predict steps fails → print and continue.
%   - No ECG channel, MTEO failure, or empty beatInfo → skip.
%   - Missing required model features → skip (prevent wrong predictions).
% Change log:
%   2025-08-30: Added unified header comments; logic unchanged.
% ========================================================================

clc; clear; close all;

% Directory setup
rootDir = pwd;
edfDir = fullfile(rootDir, 'shhs','polysomnography','edfs','shhs1');
if ~isfolder(edfDir)
    error('Directory not found: %s', edfDir);
end

% Model file
modelFile = fullfile(rootDir, 'results','trainedClassifier_latest.mat');
if ~exist(modelFile, 'file')
    error('Model file not found: %s', modelFile);
end

% Ensure paths are available
addpath(genpath(rootDir));

% Read survival mapping (nsrrid -> vital), vital: 0=Dead, 1=Alive
survivalMap = containers.Map('KeyType','char','ValueType','double');
try
    survCsv = fullfile(rootDir, 'shhs','datasets','shhs-cvd-summary-dataset-0.21.0.csv');
    if exist(survCsv, 'file')
        optsSurv = detectImportOptions(survCsv);
        Tsurv = readtable(survCsv, optsSurv);
        namesSurv = lower(Tsurv.Properties.VariableNames);
        idNsrrid = find(strcmp(namesSurv, 'nsrrid'), 1);
        idVital  = find(strcmp(namesSurv, 'vital'), 1);
        if ~isempty(idNsrrid) && ~isempty(idVital)
            nsrridCol = Tsurv.(Tsurv.Properties.VariableNames{idNsrrid});
            vitalCol  = Tsurv.(Tsurv.Properties.VariableNames{idVital});
            % Normalize to digit strings and doubles
            if iscell(nsrridCol)
                nsrridStr = cellfun(@(x) regexprep(char(x), '\\D', ''), nsrridCol, 'UniformOutput', false);
            else
                nsrridStr = regexprep(cellstr(string(nsrridCol)), '\\D', '');
            end
            if iscell(vitalCol)
                vitalNum = nan(numel(vitalCol),1);
                for kk = 1:numel(vitalCol)
                    vitalNum(kk) = str2double(string(vitalCol{kk}));
                end
            else
                vitalNum = double(vitalCol);
            end
            validMask = ~cellfun('isempty', nsrridStr) & ~isnan(vitalNum);
            for kk = 1:numel(nsrridStr)
                if validMask(kk)
                    survivalMap(nsrridStr{kk}) = vitalNum(kk);
                end
            end
        else
            fprintf('Warning: Survival CSV missing nsrrid or vital column; skipped survival mapping.\n');
        end
    else
        fprintf('Warning: Survival CSV not found: %s\n', survCsv);
    end
catch ME
    fprintf('Failed to read survival information: %s\n', ME.message);
    survivalMap = containers.Map('KeyType','char','ValueType','double');
end

% List EDFs
edfFiles = dir(fullfile(edfDir, '*.edf'));
fprintf('Found EDF files in %s: %d.\n', edfDir, numel(edfFiles));

% Fixed sampling rate
fs = 125;

% PVC probability threshold override ([] uses model's internal threshold).
% Increasing this value raises PVC precision and reduces recall; suggested 0.85~0.98.
pvcThresholdOverride = []; % Example: higher threshold to suppress false positives

% Fields required by the latest _info.mat (to decide "latest" and skip re-processing)
latestRequiredFields = {'predPVCIndices','patientVital','fs','recordNumSamples','rGlobalAll','isPVCBeat','qrs_dur_vec','r_amp_vec','sqi_vec','t_amp_vec','tGlobalIndices'};

% Load model (compatible with two save formats)
loadedModel = load(modelFile);
if isfield(loadedModel, 'trainedModelPackage')
    modelPkg = loadedModel.trainedModelPackage;
    trainedClassifier = modelPkg.trainedClassifier;
else
    trainedClassifier = loadedModel.trainedClassifier;
end

for iFile = 1:numel(edfFiles)
    edfPath = fullfile(edfFiles(iFile).folder, edfFiles(iFile).name);
    [~, recBase, ~] = fileparts(edfFiles(iFile).name);
    fprintf('\n=== [%d/%d] Processing %s ===\n', iFile, numel(edfFiles), edfFiles(iFile).name);

    % If _info.mat exists and is "latest", skip this record
    saveName = sprintf('%s_info.mat', recBase);
    savePath = fullfile(edfFiles(iFile).folder, saveName);
    [hasLatest, reasonOld] = local_is_latest_info_mat(savePath, latestRequiredFields);
    if hasLatest
        fprintf('  ✓ Latest _info.mat exists, skip: %s\n', savePath);
        continue;
    elseif exist(savePath, 'file')
        fprintf('  Outdated _info.mat found, will regenerate. Reason: %s\n', reasonOld);
    end

    % Read EDF
    try
        TT = edfread(edfPath);
    catch ME
        fprintf('  Skip: edfread failed: %s\n', ME.message);
        continue;
    end

    % Find ECG channel (prefer exact 'ECG', else fuzzy match)
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

    % Flatten to a column vector
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
        fprintf('  Skip: Unsupported ECG channel type: %s\n', class(ecgCol));
        continue;
    end
    ecg = double(ecg);
    fprintf('  ECG channel=%s, samples=%d, fs=%d Hz\n', ecgVarName, numel(ecg), fs);
    recordNumSamples = numel(ecg);

    % Filtering: skip notch (SHHS1 already notch-filtered at 60 Hz)
    try
        [ecgFiltered, ~] = ecgFilter(ecg, fs, 2, 0);
    catch ME
        fprintf('  Skip: filtering failed: %s\n', ME.message);
        continue;
    end

    % Beat detection (no external annotations)
    ATRTIMED = [];
    ANNOTD = {};
    try
        [~, beatInfo, stats] = detectAndClassifyHeartbeats(ecgFiltered, ATRTIMED, ANNOTD, fs);
    catch ME
        fprintf('  Skip: beat detection error: %s\n', ME.message);
        continue;
    end

    if isfield(stats,'mteo_failed') && stats.mteo_failed
        fprintf('  Skip: MTEO(Q/S) detection failed.\n');
        continue;
    end
    if isempty(beatInfo)
        fprintf('  Skip: no valid beats detected.\n');
        continue;
    end

    % Feature extraction (keep all beats; impute missing values before prediction)
    [featureTable, ~] = extractHeartbeatFeatures(beatInfo, fs);
    if isempty(featureTable) || height(featureTable) == 0
        fprintf('  Skip: feature extraction is empty.\n');
        continue;
    end

    % Prediction
    try
        % Select variables required by the model
        if isfield(trainedClassifier, 'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
            reqVars = trainedClassifier.RequiredVariables;
        else
            reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
        end
        missingVars = setdiff(reqVars, featureTable.Properties.VariableNames);
        if ~isempty(missingVars)
            fprintf('  Skip: features missing model required variables: %s\n', strjoin(missingVars, ', '));
            continue;
        end
        inputForPredict = featureTable(:, reqVars);

        % Median-impute missing features per column to ensure each beat is predictable
        for vv = 1:numel(reqVars)
            colName = reqVars{vv};
            colData = inputForPredict.(colName);
            if ~isnumeric(colData)
                continue;
            end
            nanMask = isnan(colData);
            if any(nanMask)
                medVal = median(colData(~nanMask));
                if isempty(medVal) || isnan(medVal)
                    medVal = 0; % fallback default
                end
                colData(nanMask) = medVal;
                inputForPredict.(colName) = colData;
            end
        end
        % Robust scalar threshold override check
        useOverride = false;
        if exist('pvcThresholdOverride','var')
            if ~isempty(pvcThresholdOverride) && isscalar(pvcThresholdOverride)
                thr = double(pvcThresholdOverride); %#ok<NASGU>
                thr = thr(1);
                if isfinite(thr) && thr >= 0 && thr <= 1
                    useOverride = true;
                    pvcThresholdOverride = thr; % enforce scalar value
                end
            end
        end
        if useOverride && isfield(trainedClassifier, 'ClassificationEnsemble')
            try
                model = trainedClassifier.ClassificationEnsemble;
                [~, rawScores] = predict(model, inputForPredict);
                classNames = model.ClassNames; % {'Other','PVC'}
                if iscell(classNames)
                    pvcIdx = find(strcmp(classNames, 'PVC'), 1);
                else
                    pvcIdx = find(classNames == 'PVC', 1);
                end
                if isempty(pvcIdx)
                    pvcIdx = min(size(rawScores,2), 2);
                end
                pvcScore = rawScores(:, pvcIdx);
                predLabels = repmat({'Other'}, size(pvcScore,1), 1);
                predLabels(pvcScore >= pvcThresholdOverride) = {'PVC'};
                fprintf('  Predicting with threshold override pvcThreshold=%.2f.\n', pvcThresholdOverride);
            catch
                % If threshold override path fails, fallback to model's built-in predictFcn
                try
                    [predLabels, ~] = trainedClassifier.predictFcn(inputForPredict);
                catch
                    predLabels = trainedClassifier.predictFcn(inputForPredict);
                end
            end
        else
            % Use model's internal threshold/rules
            try
                [predLabels, ~] = trainedClassifier.predictFcn(inputForPredict);
            catch
                predLabels = trainedClassifier.predictFcn(inputForPredict);
            end
        end
        if ~iscell(predLabels), predLabels = cellstr(predLabels); end
    catch ME
        fprintf('  Skip: prediction failed: %s\n', ME.message);
        continue;
    end

    numBeats = numel(beatInfo);
    numPVC = sum(strcmp(predLabels, 'PVC'));
    numOther = sum(strcmp(predLabels, 'Other'));
    fprintf('  Beats kept in this record: %d, predicted PVC: %d, predicted Other: %d\n', numBeats, numPVC, numOther);

    % Compute global sample indices of all R peaks and indices predicted as PVC
    rGlobal = zeros(numel(beatInfo),1);
    for kk = 1:numel(beatInfo)
        rGlobal(kk) = double(beatInfo(kk).segmentStartIndex) + double(beatInfo(kk).rIndex) - 1;
    end
    rGlobalAll = rGlobal;
    isPVCBeat = strcmp(predLabels, 'PVC');
    predPVCIndices = rGlobal(isPVCBeat);

    % Precompute per-beat quantities required by generate stage (compact)
    numBeats = numel(beatInfo);
    qrs_dur_vec = nan(numBeats,1);
    r_amp_vec   = nan(numBeats,1);
    t_amp_vec   = nan(numBeats,1);
    tGlobalIndices = nan(numBeats,1);
    % SQI vector (default true if beatInfo lacks this field)
    try
        sqi_vec = arrayfun(@(b) (isfield(b,'sqiIsGood') && ~isempty(b.sqiIsGood) && logical(b.sqiIsGood)), beatInfo);
        sqi_vec = logical(sqi_vec(:));
    catch
        sqi_vec = true(numBeats,1);
    end
    for ii = 1:numBeats
        b = beatInfo(ii);
        if isfield(b,'qrsOnIndex') && isfield(b,'qrsOffIndex') && isfinite(b.qrsOnIndex) && isfinite(b.qrsOffIndex) && b.qrsOffIndex > b.qrsOnIndex
            qrs_dur_vec(ii) = (double(b.qrsOffIndex) - double(b.qrsOnIndex))/fs;
        end
        if isfield(b,'segment') && isfield(b,'rIndex') && ~isempty(b.segment) && isfinite(b.rIndex) && b.rIndex>=1 && b.rIndex<=numel(b.segment)
            r_amp_vec(ii) = double(b.segment(b.rIndex));
        end
        if isfield(b,'tIndex') && isfinite(b.tIndex) && b.tIndex>=1 && isfield(b,'segment') && ~isempty(b.segment) && b.tIndex<=numel(b.segment)
            t_amp_vec(ii) = abs(double(b.segment(b.tIndex)));
            if isfield(b,'segmentStartIndex')
                tGlobalIndices(ii) = double(b.segmentStartIndex) + double(b.tIndex) - 1;
            end
        end
    end

    % Get patient vital status
    patientVital = NaN;
    try
        idStr = regexprep(extractAfter(recBase, 'shhs1-'), '\\D', '');
        if isKey(survivalMap, idStr)
            patientVital = survivalMap(idStr);
        end
    catch
        % Keep as NaN
    end

    % Save to the same EDF directory (compact set required by generate)
    saveName = sprintf('%s_info.mat', recBase);
    savePath = fullfile(edfFiles(iFile).folder, saveName);
    % Metadata: used to determine "latest version"
    infoSchemaVersion = 3;              % Current _info.mat schema version
    infoGenerator = 'predict_shhs1_all.m';
    infoGeneratedAt = char(datetime("now","Format","yyyy-MM-dd HH:mm:ss"));
    try
        save(savePath, 'predPVCIndices', 'patientVital', 'fs', 'recordNumSamples', 'rGlobalAll', 'isPVCBeat', 'qrs_dur_vec', 'r_amp_vec', 'sqi_vec', 't_amp_vec', 'tGlobalIndices', 'infoSchemaVersion', 'infoGenerator', 'infoGeneratedAt');
        fprintf('  ✓ Saved: %s\n', savePath);
    catch ME
        fprintf('  Save failed: %s\n', ME.message);
    end
end

fprintf('\nAll processing done.\n');

% ===================== Local functions in this script =====================
function [isLatest, reason] = local_is_latest_info_mat(matPath, requiredFields)
% Determine whether the specified _info.mat contains the latest required fields, with basic consistency checks
    isLatest = false;
    reason = '';
    if exist(matPath, 'file') ~= 2
        reason = 'File not found';
        return;
    end
    try
        finfo = whos('-file', matPath);
        vars = {finfo.name};
        miss = setdiff(requiredFields, vars);
        if ~isempty(miss)
            reason = ['Missing fields: ' strjoin(miss, ', ')];
            return;
        end
        % Further minimal consistency checks (types/lengths)
        S = load(matPath, 'fs','recordNumSamples','rGlobalAll','isPVCBeat','qrs_dur_vec','r_amp_vec','sqi_vec','t_amp_vec','tGlobalIndices','patientVital','predPVCIndices','infoSchemaVersion');
        if ~isfield(S,'fs') || ~isnumeric(S.fs) || ~isscalar(S.fs) || ~isfinite(S.fs)
            reason = 'fs is not a numeric scalar'; return;
        end
        if ~isfield(S,'recordNumSamples') || ~isnumeric(S.recordNumSamples) || ~isscalar(S.recordNumSamples) || ~isfinite(S.recordNumSamples)
            reason = 'recordNumSamples is not a numeric scalar'; return;
        end
        L = numel(S.rGlobalAll);
        if any([numel(S.isPVCBeat), numel(S.qrs_dur_vec), numel(S.r_amp_vec), numel(S.sqi_vec), numel(S.t_amp_vec), numel(S.tGlobalIndices)] ~= L)
            reason = 'Per-beat vector lengths are inconsistent'; return;
        end
        % If a version number exists and is older, treat as not latest
        if isfield(S,'infoSchemaVersion')
            try
                ver = double(S.infoSchemaVersion);
            catch
                ver = NaN;
            end
            if isfinite(ver) && ver < 3
                reason = sprintf('infoSchemaVersion=%g older than 3', ver);
                return;
            end
        end
        isLatest = true;
    catch ME
        reason = sprintf('Exception during reading/validation: %s', ME.message);
        isLatest = false;
    end
end
