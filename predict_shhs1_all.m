% File: predict_shhs1_all.m
% Type: Script
% Description:
%   Iterate all SHHS1 EDFs, segment beats using existing pipeline, predict beat types using the trained model,
%   and save predictions and features to the corresponding directories.
% Usage:
%   Run this script directly (requires results/trainedClassifier_latest.mat).
% Inputs/Dependencies:
%   - Data dir: shhs/polysomnography/edfs/shhs1
%   - Functions: ecgFilter, detectAndClassifyHeartbeats, extractHeartbeatFeatures
% Outputs:
%   - *_beats_features_pred.mat per record
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

clc; clear; close all;

% Directory
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

% Ensure path
addpath(genpath(rootDir));

% List EDFs
edfFiles = dir(fullfile(edfDir, '*.edf'));
fprintf('Found %d EDF files in %s.\n', numel(edfFiles), edfDir);

% Fixed sampling rate
fs = 125;

% PVC probability threshold override ([] means use model default).
% Higher value increases precision but reduces recall; suggested 0.85~0.98.
pvcThresholdOverride = 0.92; % example: higher threshold to reduce FPs

% Load model (either packaging styles)
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

    % Read EDF
    try
        TT = edfread(edfPath);
    catch ME
        fprintf('  Skip: edfread failed: %s\n', ME.message);
        continue;
    end

    % Find ECG channel ('ECG' exact else fuzzy)
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

    % Flatten to column
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
    fprintf('  ECG=%s, samples=%d, fs=%d Hz\n', ecgVarName, numel(ecg), fs);

    % Filtering: skip notch (SHHS1 already notched at 60 Hz)
    try
        [ecgFiltered, ~] = ecgFilter(ecg, fs, 2, 0);
    catch ME
        fprintf('  Skip: filtering failed: %s\n', ME.message);
        continue;
    end

    % Beat detection (no annotations)
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
        fprintf('  Skip: no valid beats.\n');
        continue;
    end

    % Feature extraction (keep all beats; impute later for prediction)
    [featureTable, ~] = extractHeartbeatFeatures(beatInfo, fs);
    if isempty(featureTable) || height(featureTable) == 0
        fprintf('  Skip: feature extraction is empty.\n');
        continue;
    end

    % Prediction
    try
        % Select required variables by model
        if isfield(trainedClassifier, 'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
            reqVars = trainedClassifier.RequiredVariables;
        else
            reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
        end
        missingVars = setdiff(reqVars, featureTable.Properties.VariableNames);
        if ~isempty(missingVars)
            fprintf('  Skip: features missing required variables: %s\n', strjoin(missingVars, ', '));
            continue;
        end
        inputForPredict = featureTable(:, reqVars);

        % Median impute missing numeric features to ensure predictability
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
                    medVal = 0; % fallback
                end
                colData(nanMask) = medVal;
                inputForPredict.(colName) = colData;
            end
        end
        useOverride = exist('pvcThresholdOverride','var') && ~isempty(pvcThresholdOverride) && isfinite(pvcThresholdOverride) && pvcThresholdOverride >= 0 && pvcThresholdOverride <= 1;
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
                fprintf('  Using threshold override pvcThreshold=%.2f.\n', pvcThresholdOverride);
            catch
                % Fallback to model predictFcn if override fails
                try
                    [predLabels, ~] = trainedClassifier.predictFcn(inputForPredict);
                catch
                    predLabels = trainedClassifier.predictFcn(inputForPredict);
                end
            end
        else
            % Use model default threshold/logic
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

    % Summarize and sync to allBeatInfo (cover all beats)
    allBeatInfo = beatInfo; %#ok<NASGU>
    for k = 1:numel(allBeatInfo)
        bt = predLabels{k};
        if ~ischar(bt), bt = char(string(bt)); end
        allBeatInfo(k).beatType = bt;
        allBeatInfo(k).fs = fs;
        allBeatInfo(k).originalRecordName = recBase;
        allBeatInfo(k).originalDatabaseName = 'SHHS1';
    end

    numBeats = numel(allBeatInfo);
    numPVC = sum(strcmp(predLabels, 'PVC'));
    numOther = sum(strcmp(predLabels, 'Other'));
    fprintf('  Kept beats: %d, predicted PVC: %d, Other: %d\n', numBeats, numPVC, numOther);

    % Save next to EDF
    saveName = sprintf('%s_beats_features_pred.mat', recBase);
    savePath = fullfile(edfFiles(iFile).folder, saveName);
    try
        edfPathSaved = edfPath; %#ok<NASGU>
        ecgVarNameSaved = ecgVarName; %#ok<NASGU>
        statsSaved = stats; %#ok<NASGU>
        predLabelsSaved = predLabels; %#ok<NASGU>
        featureTable_raw = featureTable; %#ok<NASGU>
        featureTable_imputed = inputForPredict; %#ok<NASGU>
        save(savePath, 'allBeatInfo', 'featureTable_raw', 'featureTable_imputed', 'predLabelsSaved', 'statsSaved', 'fs', 'ecgVarNameSaved', 'edfPathSaved');
        fprintf('  ✓ Saved: %s\n', savePath);
    catch ME
        fprintf('  Save failed: %s\n', ME.message);
    end
end

fprintf('\nAll done.\n');


