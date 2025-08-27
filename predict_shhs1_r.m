% File: predict_shhs1_r.m
% Type: Script
% Description:
%   Process only SHHS1 records with rpoints annotations, predict beat types and evaluate by
%   nearest-neighbor matching to annotations, and save *_info.mat (with predPVCIndices and patientVital).
% Usage:
%   Run directly (requires results/trainedClassifier_latest.mat and rpoints CSV present).
% Inputs/Dependencies:
%   - Data dirs: shhs/polysomnography/edfs/shhs1, shhs/polysomnography/annotations-rpoints/shhs1
%   - Functions: ecgFilter, detectAndClassifyHeartbeats, extractHeartbeatFeatures
% Outputs:
%   - *_beats_features_pred.mat and *_info.mat per record
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26


clc; clear; close all;

% Directory settings
rootDir = pwd;
edfDir = fullfile(rootDir, 'shhs','polysomnography','edfs','shhs1');
annDir = fullfile(rootDir, 'shhs','polysomnography','annotations-rpoints','shhs1');
if ~isfolder(edfDir)
    error('EDF directory not found: %s', edfDir);
end
if ~isfolder(annDir)
    error('rpoints directory not found: %s', annDir);
end

% Model file
modelFile = fullfile(rootDir, 'results','trainedClassifier_latest.mat');
if ~exist(modelFile, 'file')
    error('Model file not found: %s', modelFile);
end

% Ensure path is available
addpath(genpath(rootDir));

% Read survival info map (nsrrid -> vital), vital: 0=Dead, 1=Alive
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
            % Normalize to numeric strings and double
            if iscell(nsrridCol)
                nsrridStr = cellfun(@(x) regexprep(char(x), '\D', ''), nsrridCol, 'UniformOutput', false);
            else
                nsrridStr = regexprep(cellstr(string(nsrridCol)), '\D', '');
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
            fprintf('Warning: survival CSV missing nsrrid or vital; skip survival map.\n');
        end
    else
        fprintf('Warning: survival CSV not found: %s\n', survCsv);
    end
catch ME
    fprintf('Failed to read survival info: %s\n', ME.message);
    survivalMap = containers.Map('KeyType','char','ValueType','double');
end

% List records with rpoints annotations (shhs1-XXXXX-rpoint.csv → shhs1-XXXXX)
annFiles = dir(fullfile(annDir, '*-rpoint.csv'));
if isempty(annFiles)
    error('rpoints directory is empty: %s', annDir);
end

recordBases = cell(numel(annFiles), 1);
for i = 1:numel(annFiles)
    nm = annFiles(i).name;
    % Remove suffix "-rpoint.csv"
    if endsWith(nm, '-rpoint.csv', 'IgnoreCase', true)
        recordBases{i} = extractBefore(nm, strlength(nm) - strlength('-rpoint.csv') + 1);
    else
        recordBases{i} = ''; % skip non-matching name
    end
end
recordBases = recordBases(~cellfun('isempty', recordBases));
recordBases = unique(recordBases);

% Map to existing EDFs
edfPaths = {};
for i = 1:numel(recordBases)
    edfPathCand = fullfile(edfDir, [recordBases{i} '.edf']);
    if exist(edfPathCand, 'file')
        edfPaths{end+1,1} = edfPathCand; %#ok<AGROW>
    end
end

fprintf('Found %d EDF records with rpoints in %s (total %d rpoints files).\n', numel(edfPaths), edfDir, numel(annFiles));
if isempty(edfPaths)
    fprintf('No records to process. Exiting.\n');
    return;
end

% Fixed sampling rate
fs = 125;

% PVC probability threshold override ([] uses model default). Higher value → more precision, less recall.
pvcThresholdOverride = 0.92;

% Load model (compatible with both saving styles)
loadedModel = load(modelFile);
if isfield(loadedModel, 'trainedModelPackage')
    modelPkg = loadedModel.trainedModelPackage; %#ok<NASGU>
    trainedClassifier = modelPkg.trainedClassifier;
else
    trainedClassifier = loadedModel.trainedClassifier;
end

for iFile = 1:numel(edfPaths)
    edfPath = edfPaths{iFile};
    [edfFolder, edfName, ~] = fileparts(edfPath);
    recBase = edfName; % e.g., shhs1-XXXXXX
    fprintf('\n=== [%d/%d] Processing %s ===\n', iFile, numel(edfPaths), [edfName '.edf']);

    % Before prediction: print PVC/Other stats from annotation CSV (Type: 0=Artifact,1=NSR,2=VE,3=SVE)
    try
        annCsv = fullfile(annDir, [recBase '-rpoint.csv']);
        if exist(annCsv, 'file')
            optsAnn = detectImportOptions(annCsv);
            Tann = readtable(annCsv, optsAnn);
            namesAnn = lower(Tann.Properties.VariableNames);
            idTypeAnn = find(strcmp(namesAnn, 'type'), 1);
            if ~isempty(idTypeAnn)
                typeRaw = Tann.(Tann.Properties.VariableNames{idTypeAnn});
                if iscell(typeRaw)
                    typeVec = nan(size(typeRaw));
                    for kk = 1:numel(typeRaw), typeVec(kk) = str2double(string(typeRaw{kk})); end
                else
                    typeVec = double(typeRaw);
                end
                mValid = isfinite(typeVec) & (typeVec ~= 0);
                typeValid = typeVec(mValid);
                numPVC_ann = sum(typeValid == 2);
                numOther_ann = sum((typeValid == 1) | (typeValid == 3));
                fprintf('  Annotation stats: PVC=%d, Other=%d\n', numPVC_ann, numOther_ann);
            else
                fprintf('  Skip stats: rpoints missing Type field.\n');
            end
        else
            fprintf('  Skip stats: rpoints file not found.\n');
        end
    catch ME
        fprintf('  Stats failed: %s\n', ME.message);
    end

    % Read EDF (ECG only)
    try
        TT = edfread(edfPath);
    catch ME
        fprintf('  Skip: edfread failed: %s\n', ME.message);
        continue;
    end

    % Find ECG channel (exact 'ECG' else fuzzy)
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

    % Flatten to column vector
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
        fprintf('  Skip: unsupported ECG channel type: %s\n', class(ecgCol));
        continue;
    end
    ecg = double(ecg);
    fprintf('  ECG=%s, samples=%d, fs=%d Hz\n', ecgVarName, numel(ecg), fs);

    % Note: In some leads the entire signal may be inverted; optionally invert if needed
    % ecg = -ecg;

    % Filtering: skip notch (SHHS1 already notched at 60 Hz)
    try
        [ecgFiltered, ~] = ecgFilter(ecg, fs, 2, 0);
    catch ME
        fprintf('  Skip: filtering failed: %s\n', ME.message);
        continue;
    end

    % Beat detection (no annotations; rpoints only for record selection)
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
        fprintf('  Skip: feature extraction empty.\n');
        continue;
    end

    % Prediction
    try
        % Select model required variables
        if isfield(trainedClassifier, 'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
            reqVars = trainedClassifier.RequiredVariables;
        else
            reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
        end
        missingVars = setdiff(reqVars, featureTable.Properties.VariableNames);
        if ~isempty(missingVars)
            fprintf('  Skip: features missing required vars: %s\n', strjoin(missingVars, ', '));
            continue;
        end
        inputForPredict = featureTable(:, reqVars);

        % Median-impute missing numeric features
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

    % Summarize and sync to allBeatInfo (all beats)
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

    % Compute global sample indices of predicted PVC R-peaks
    rGlobal = zeros(numel(allBeatInfo),1);
    for kk = 1:numel(allBeatInfo)
        rGlobal(kk) = double(allBeatInfo(kk).segmentStartIndex) + double(allBeatInfo(kk).rIndex) - 1;
    end
    predPVCIndices = rGlobal(strcmp(predLabels, 'PVC'));

    % === Per-record evaluation based on rpoints (nearest matching) ===
    try
        annCsv = fullfile(annDir, [recBase '-rpoint.csv']);
        if exist(annCsv, 'file')
            optsEval = detectImportOptions(annCsv);
            TannEval = readtable(annCsv, optsEval);
            namesEval = lower(TannEval.Properties.VariableNames);
            idType = find(strcmp(namesEval, 'type'), 1);
            idSec  = find(strcmp(namesEval, 'seconds'), 1);
            if ~isempty(idType) && ~isempty(idSec)
                typeRawEval = TannEval.(TannEval.Properties.VariableNames{idType});
                secRawEval  = TannEval.(TannEval.Properties.VariableNames{idSec});
                % Convert columns
                if iscell(typeRawEval)
                    typeVecEval = nan(size(typeRawEval));
                    for kk = 1:numel(typeRawEval), typeVecEval(kk) = str2double(string(typeRawEval{kk})); end
                else
                    typeVecEval = double(typeRawEval);
                end
                if iscell(secRawEval)
                    secVecEval = nan(size(secRawEval));
                    for kk = 1:numel(secRawEval), secVecEval(kk) = str2double(string(secRawEval{kk})); end
                else
                    secVecEval = double(secRawEval);
                end
                % Filter valid ground-truth (exclude Type==0)
                mValid = isfinite(secVecEval) & ~isnan(secVecEval) & ~isnan(typeVecEval) & (typeVecEval ~= 0);
                secGT = secVecEval(mValid);
                typeGTnum = typeVecEval(mValid);
                annGT = cell(numel(typeGTnum),1);
                for kk = 1:numel(typeGTnum)
                    if typeGTnum(kk) == 2
                        annGT{kk} = 'PVC';
                    else
                        annGT{kk} = 'Other';
                    end
                end
                % Convert predicted beats' global R indices to seconds (beatInfo has segmentStartIndex and rIndex)
                rGlobal = zeros(numel(allBeatInfo),1);
                for kk = 1:numel(allBeatInfo)
                    rGlobal(kk) = double(allBeatInfo(kk).segmentStartIndex) + double(allBeatInfo(kk).rIndex) - 1;
                end
                rSec = rGlobal(:) / fs;
                % Nearest matching tolerance: ±0.1 s
                tolSec = 0.10;
                yTrue = {};
                yPred = {};
                for kk = 1:numel(secGT)
                    [minDiff, idxP] = min(abs(rSec - secGT(kk)));
                    if ~isempty(idxP) && ~isnan(minDiff) && minDiff <= tolSec
                        yTrue{end+1,1} = annGT{kk}; %#ok<AGROW>
                        yPred{end+1,1} = predLabels{idxP}; %#ok<AGROW>
                    end
                end
                if ~isempty(yTrue)
                    % Text confusion matrix and metrics
                    classes = unique(yTrue);
                    nC = numel(classes);
                    C = zeros(nC, nC);
                    for ii = 1:nC
                        for jj = 1:nC
                            C(ii,jj) = sum(strcmp(yTrue, classes{ii}) & strcmp(yPred, classes{jj}));
                        end
                    end
                    fprintf('  — Text Confusion Matrix (record %s) —\n', recBase);
                    hdr = sprintf('%-12s', 'Actual\\Pred');
                    for ii = 1:nC, hdr = [hdr, sprintf('%12s', classes{ii})]; end %#ok<AGROW>
                    fprintf('%s\n', hdr);
                    fprintf('%s\n', repmat('-', 1, length(hdr)));
                    for ii = 1:nC
                        row = sprintf('%-12s', classes{ii});
                        for jj = 1:nC
                            row = [row, sprintf('%12d', C(ii,jj))]; %#ok<AGROW>
                        end
                        fprintf('%s\n', row);
                    end
                    % Metrics
                    for ii = 1:nC
                        TP = C(ii,ii);
                        FP = sum(C(:,ii)) - TP;
                        FN = sum(C(ii,:)) - TP;
                        TN = sum(C(:)) - TP - FP - FN;
                        if (TP + FN) > 0, sens = TP/(TP+FN)*100; else, sens = 0; end
                        if (TN + FP) > 0, spec = TN/(TN+FP)*100; else, spec = 0; end
                        if (TP + FP) > 0, prec = TP/(TP+FP)*100; else, prec = 0; end
                        if (prec + sens) > 0, f1 = 2*(prec*sens)/(prec+sens); else, f1 = 0; end
                        fprintf('  Class %s: TPR=%.2f%%, FNR=%.2f%%, Specificity=%.2f%%, Precision=%.2f%%, F1=%.2f\n', ...
                            classes{ii}, sens, 100-sens, spec, prec, f1);
                    end
                else
                    fprintf('  Eval note: no matches within ±0.1 s.\n');
                end
            else
                fprintf('  Eval skipped: rpoints missing Type/seconds fields.\n');
            end
        else
            fprintf('  Eval skipped: rpoints file not found.\n');
        end
    catch ME
        fprintf('  Eval failed: %s\n', ME.message);
    end

    % Retrieve patient survival info
    patientVital = NaN;
    try
        idStr = regexprep(extractAfter(recBase, 'shhs1-'), '\D', '');
        if isKey(survivalMap, idStr)
            patientVital = survivalMap(idStr);
        end
    catch
        % keep NaN
    end

    % Save next to EDF (predPVCIndices & patientVital only)
    saveName = sprintf('%s_info.mat', recBase);
    savePath = fullfile(edfFolder, saveName);
    try
        save(savePath, 'predPVCIndices', 'patientVital');
        fprintf('  ✓ Saved: %s\n', savePath);
    catch ME
        fprintf('  Save failed: %s\n', ME.message);
    end
end

fprintf('\nAll done.\n');


