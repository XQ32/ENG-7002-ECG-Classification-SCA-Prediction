% File: generate_beats_classifier.m
% Type: Script
% Description:
%   Build beat feature tables from SHHS1 EDF + rpoints, stratify train/test,
%   train PVC/Other classifier and save model and feature data.
% Usage:
%   Run this script directly (requires SHHS directory structure and dependencies on path).
% Inputs/Dependencies:
%   - Data dirs: shhs/polysomnography/edfs/shhs1,
%                shhs/polysomnography/annotations-rpoints/shhs1
%   - Functions: ecgFilter, detectAndClassifyHeartbeats,
%                extractHeartbeatFeatures, trainBeatsClassifier
% Outputs:
%   - results/trainingFeatureTable.mat, results/testingFeatureTable.mat
%   - results/trainedClassifier_latest.mat (and timestamped copy)
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26


clc; clear; close all;

try
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

    % Path setup
    addpath(genpath(rootDir));

    % Split/train config
    splitOptions = struct();
    splitOptions.test_split_ratio = 0.2;  % test ratio
    splitOptions.random_seed = 42;        % random seed
    splitOptions.stratified = true;       % stratified sampling

    % Process only first N annotated EDFs (<=0 means all)
    processFirstN = 100;

    % Training options
    model_options = struct();
    model_options.enableCostMatrix = true;
    % High precision: increase Other->PVC cost, reduce PVC->Other cost
    model_options.costFP = 5;  % false positive cost (Other -> PVC)
    model_options.costFN = 1;  % false negative cost (PVC -> Other)
    model_options.enableThresholdMoving = true;
    model_options.pvcThreshold = 0.95; % higher threshold to reduce FPs

    % List records with rpoints annotations (shhs1-XXXXX-rpoint.csv → shhs1-XXXXX)
    annFiles = dir(fullfile(annDir, '*-rpoint.csv'));
    if isempty(annFiles)
        error('rpoints directory is empty: %s', annDir);
    end

    recordBases = cell(numel(annFiles), 1);
    for i = 1:numel(annFiles)
        nm = annFiles(i).name;
        if endsWith(nm, '-rpoint.csv', 'IgnoreCase', true)
            recordBases{i} = extractBefore(nm, strlength(nm) - strlength('-rpoint.csv') + 1);
        else
            recordBases{i} = '';
        end
    end
    recordBases = recordBases(~cellfun('isempty', recordBases));
    recordBases = unique(recordBases);

    % Filter to existing EDFs
    edfPaths = cell(0,1);
    annPaths = cell(0,1);
    for i = 1:numel(recordBases)
        edfPathCand = fullfile(edfDir, [recordBases{i} '.edf']);
        annPathCand = fullfile(annDir, [recordBases{i} '-rpoint.csv']);
        if exist(edfPathCand, 'file') && exist(annPathCand, 'file')
            edfPaths{end+1,1} = edfPathCand; %#ok<AGROW>
            annPaths{end+1,1} = annPathCand; %#ok<AGROW>
        end
    end

    % Truncate if only first N are processed
    if exist('processFirstN','var') && isnumeric(processFirstN) && ~isempty(processFirstN) && processFirstN > 0
        maxN = min(processFirstN, numel(edfPaths));
        edfPaths = edfPaths(1:maxN);
        annPaths = annPaths(1:maxN);
        fprintf('Found %d EDF records with rpoints in %s (total %d rpoints files). Processing first %d.\n', numel(recordBases), edfDir, numel(annFiles), maxN);
    else
        fprintf('Found %d EDF records with rpoints in %s (total %d rpoints files).\n', numel(edfPaths), edfDir, numel(annFiles));
    end
    if isempty(edfPaths)
        fprintf('No records to process. Exiting.\n');
        return;
    end

    % Keep only feature tables (per-record temporary storage)
    features_by_record = cell(numel(edfPaths), 1);
    numKept = 0;

    for iFile = 1:numel(edfPaths)
        edfPath = edfPaths{iFile};
        annPath = annPaths{iFile};
        [edfFolder, edfName, ~] = fileparts(edfPath);
        recBase = edfName; % e.g., shhs1-XXXXXX
        fprintf('\n=== [%d/%d] Processing %s ===\n', iFile, numel(edfPaths), [edfName '.edf']);

        % Read EDF (ECG only)
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

        % Read rpoints CSV: seconds/Type/samplingrate
        try
            opts = detectImportOptions(annPath);
            Tann = readtable(annPath, opts);
        catch ME
            fprintf('  Skip: reading rpoints failed: %s\n', ME.message);
            continue;
        end
        vnames = Tann.Properties.VariableNames;
        vnamesLower = lower(vnames);
        idxType = find(strcmp(vnamesLower, 'type'), 1);
        idxSec  = find(strcmp(vnamesLower, 'seconds'), 1);
        idxSR   = find(strcmp(vnamesLower, 'samplingrate'), 1);
        if isempty(idxType) || isempty(idxSec)
            fprintf('  Skip: rpoints missing required fields (Type/seconds).\n');
            continue;
        end

        % Sampling rate (prefer CSV; fallback 125 Hz)
        fsLocal = 125;
        if ~isempty(idxSR)
            srRaw = Tann.(vnames{idxSR});
            if isnumeric(srRaw)
                srVec = double(srRaw(:));
            elseif iscell(srRaw)
                srVec = nan(numel(srRaw),1);
                for ii = 1:numel(srRaw), srVec(ii) = str2double(string(srRaw{ii})); end
            else
                srVec = double(srRaw);
            end
            srVec = srVec(isfinite(srVec) & ~isnan(srVec));
            if ~isempty(srVec)
                fsCand = round(median(srVec,'omitnan'));
                if isfinite(fsCand) && ~isnan(fsCand) && fsCand > 0
                    fsLocal = fsCand;
                end
            end
        end
        fprintf('  ECG channel=%s, samples=%d, fs=%d Hz (from rpoints)\n', ecgVarName, numel(ecg), fsLocal);

        % seconds/Type conversion
        typeRaw = Tann.(vnames{idxType});
        secRaw  = Tann.(vnames{idxSec});
        if iscell(typeRaw)
            typeVec = nan(size(typeRaw));
            for ii = 1:numel(typeRaw), typeVec(ii) = str2double(string(typeRaw{ii})); end
        else
            typeVec = double(typeRaw);
        end
        if iscell(secRaw)
            secVec = nan(size(secRaw));
            for ii = 1:numel(secRaw), secVec(ii) = str2double(string(secRaw{ii})); end
        else
            secVec = double(secRaw);
        end
        % Filter: remove Type==0 (Artifact) and invalid seconds
        mask_valid = ~isnan(secVec) & isfinite(secVec) & ~isnan(typeVec) & (typeVec ~= 0);
        ATRTIMED = secVec(mask_valid);
        type_kept = typeVec(mask_valid);
        if isempty(ATRTIMED)
            fprintf('  Skip: rpoints valid annotations is empty.\n');
            continue;
        end
        % Type mapping: 2->PVC; 1/3->Other; others -> Other
        ANNOTD = cell(numel(type_kept), 1);
        for ii = 1:numel(type_kept)
            if type_kept(ii) == 2
                ANNOTD{ii} = 'PVC';
            else
                ANNOTD{ii} = 'Other';
            end
        end

        % Print annotation counts
        pvcCount = sum(strcmp(ANNOTD, 'PVC'));
        otherCount = sum(strcmp(ANNOTD, 'Other'));
        fprintf('  Annotation: PVC=%d, Other=%d (total=%d)\n', pvcCount, otherCount, pvcCount + otherCount);

        % Optional: invert entire ECG
        % ecg = -ecg;

        % Filtering: skip notch (60 Hz handled)
        try
            [ecgFiltered, ~] = ecgFilter(ecg, fsLocal, 2, 0);
        catch ME
            fprintf('  Skip: filtering failed: %s\n', ME.message);
            continue;
        end

        % Beat detection and typing (nearest R-annotation)
        try
            [~, beatInfo, stats] = detectAndClassifyHeartbeats(ecgFiltered, ATRTIMED, ANNOTD, fsLocal);
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

        % Feature extraction (keep features; drop missing)
        [featureTable_rec, ~] = extractHeartbeatFeatures(beatInfo, fsLocal);
        if isempty(featureTable_rec) || height(featureTable_rec) == 0
            fprintf('  Skip: feature extraction empty.\n');
            continue;
        end
        completeRows = ~any(ismissing(featureTable_rec), 2);
        featureTable_rec = featureTable_rec(completeRows, :);
        fprintf('  Valid beats (feature rows): %d\n', height(featureTable_rec));

        numKept = numKept + 1;
        features_by_record{numKept,1} = featureTable_rec;

        % Cleanup temporary variables for this record
        clear TT ecg ecgCol ecgVarName ecgFiltered Tann ATRTIMED ANNOTD featureTable_rec completeRows ...
              typeRaw secRaw typeVec secVec mask_valid type_kept idxType idxSec idxSR srRaw srVec fsLocal;
    end

    % Merge all features
    if numKept == 0
        fprintf('\nNo valid features extracted from any record. Exit.\n');
        return;
    end
    features_by_record = features_by_record(1:numKept);
    valid_tables = features_by_record(~cellfun('isempty', features_by_record));
    allFeatureTable = vertcat(valid_tables{:});
    fprintf('\nMerged all features. Total beats: %d\n', height(allFeatureTable));

    % Stratified random split
    fprintf('\n=== Stratified split ===\n');
    fprintf('Test ratio: %.1f%%, seed: %d\n', splitOptions.test_split_ratio*100, splitOptions.random_seed);
    rng(splitOptions.random_seed);
    all_beat_types = allFeatureTable.BeatType;
    unique_types = unique(all_beat_types);
    train_indices = [];
    test_indices = [];
    for type_idx = 1:length(unique_types)
        current_type = unique_types{type_idx};
        type_indices = find(strcmp(all_beat_types, current_type));
        type_count = length(type_indices);
        test_count = round(type_count * splitOptions.test_split_ratio);
        train_count = type_count - test_count; %#ok<NASGU>
        shuffled_indices = type_indices(randperm(type_count));
        test_indices = [test_indices; shuffled_indices(1:test_count)];
        train_indices = [train_indices; shuffled_indices(test_count+1:end)];
        fprintf('Class %s: total=%d, train=%d, test=%d (%.1f%%)\n', ...
            current_type, type_count, train_count, test_count, test_count/type_count*100);
    end
    trainingFeatureTable = allFeatureTable(train_indices, :);
    testingFeatureTable = allFeatureTable(test_indices, :);
    fprintf('\nSplit done: train=%d, test=%d\n', height(trainingFeatureTable), height(testingFeatureTable));

    % Save feature tables
    if ~exist('results','dir'), mkdir('results'); end
    try
        save(fullfile('results','trainingFeatureTable.mat'), 'trainingFeatureTable');
        fprintf('✓ Saved training feature table\n');
    catch ME
        fprintf('⚠ Failed to save training feature table: %s\n', ME.message);
    end
    try
        save(fullfile('results','testingFeatureTable.mat'), 'testingFeatureTable');
        fprintf('✓ Saved testing feature table\n');
    catch ME
        fprintf('⚠ Failed to save testing feature table: %s\n', ME.message);
    end

    % Train classifier
    if height(trainingFeatureTable) > 0
        fprintf('\n--- Train classifier ---\n');
        [trainedClassifier, validationAccuracy] = trainBeatsClassifier(trainingFeatureTable, model_options);
        fprintf('Classifier trained. CV accuracy: %.2f%%\n', validationAccuracy * 100);

        % Evaluation on test set
        if height(testingFeatureTable) > 0
            fprintf('\n--- Evaluate on test set ---\n');
            try
                [predicted_labels, prediction_scores] = trainedClassifier.predictFcn(testingFeatureTable); %#ok<NASGU>
            catch
                predicted_labels = trainedClassifier.predictFcn(testingFeatureTable);
            end
            actual_labels = testingFeatureTable.BeatType;
            if iscell(predicted_labels) && iscell(actual_labels)
                accuracy = sum(strcmp(predicted_labels, actual_labels)) / length(actual_labels) * 100;
            elseif iscell(predicted_labels) && ~iscell(actual_labels)
                accuracy = sum(strcmp(predicted_labels, cellstr(actual_labels))) / length(actual_labels) * 100;
            elseif ~iscell(predicted_labels) && iscell(actual_labels)
                accuracy = sum(strcmp(cellstr(predicted_labels), actual_labels)) / length(actual_labels) * 100;
            else
                accuracy = sum(predicted_labels == actual_labels) / length(actual_labels) * 100;
            end
            fprintf('Test accuracy: %.2f%%\n', accuracy);

            % === Confusion Matrix (figure + text) ===
            % Normalize to cellstr
            if ~iscell(actual_labels), actual_labels = cellstr(actual_labels); end
            if ~iscell(predicted_labels), predicted_labels = cellstr(predicted_labels); end
            class_names = unique(actual_labels);
            num_classes = numel(class_names);
            confusion_mat = zeros(num_classes, num_classes);
            for ii = 1:num_classes
                for jj = 1:num_classes
                    confusion_mat(ii, jj) = sum(strcmp(actual_labels, class_names{ii}) & strcmp(predicted_labels, class_names{jj}));
                end
            end

            % Plot confusion matrix (row/column summaries)
            figure('Name', 'Confusion Matrix (Test)', 'NumberTitle', 'off');
            confusionchart(confusion_mat, class_names, 'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
            title('ConfusionMatrix (Test Set)');

            % Text (absolute counts)
            fprintf('\n=== Confusion Matrix (Counts) ===\n');
            header = sprintf('%-12s', 'Actual\\Pred');
            for ii = 1:num_classes
                header = [header, sprintf('%12s', class_names{ii})]; %#ok<AGROW>
            end
            fprintf('%s\n', header);
            fprintf('%s\n', repmat('-', 1, length(header)));
            for ii = 1:num_classes
                row_str = sprintf('%-12s', class_names{ii});
                for jj = 1:num_classes
                    row_str = [row_str, sprintf('%12d', confusion_mat(ii, jj))]; %#ok<AGROW>
                end
                fprintf('%s\n', row_str);
            end

            % Text (row percentage)
            fprintf('\n=== Confusion Matrix (Row %) ===\n');
            row_sums = sum(confusion_mat, 2);
            header = sprintf('%-12s', 'Actual\\Pred');
            for ii = 1:num_classes
                header = [header, sprintf('%12s', class_names{ii})]; %#ok<AGROW>
            end
            fprintf('%s\n', header);
            fprintf('%s\n', repmat('-', 1, length(header)));
            for ii = 1:num_classes
                row_str = sprintf('%-12s', class_names{ii});
                for jj = 1:num_classes
                    if row_sums(ii) > 0
                        percentage = confusion_mat(ii, jj) / row_sums(ii) * 100;
                        row_str = [row_str, sprintf('%11.1f%%', percentage)]; %#ok<AGROW>
                    else
                        row_str = [row_str, sprintf('%11s', 'N/A')]; %#ok<AGROW>
                    end
                end
                fprintf('%s\n', row_str);
            end

            % Per-class metrics
            fprintf('\n=== Per-class Metrics ===\n');
            for ii = 1:num_classes
                TP = confusion_mat(ii, ii);
                FP = sum(confusion_mat(:, ii)) - TP;
                FN = sum(confusion_mat(ii, :)) - TP;
                TN = sum(confusion_mat(:)) - TP - FP - FN;
                if (TP + FN) > 0
                    sensitivity = TP / (TP + FN) * 100; % TPR
                else
                    sensitivity = 0;
                end
                if (TN + FP) > 0
                    specificity = TN / (TN + FP) * 100;
                else
                    specificity = 0;
                end
                if (TP + FP) > 0
                    precision = TP / (TP + FP) * 100;
                else
                    precision = 0;
                end
                if (precision + sensitivity) > 0
                    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity);
                else
                    f1_score = 0;
                end
                fprintf('Class %s: TPR=%.2f%%, FNR=%.2f%%, Specificity=%.2f%%, Precision=%.2f%%, F1=%.2f\n', ...
                    class_names{ii}, sensitivity, 100 - sensitivity, specificity, precision, f1_score);
            end
        end

        % Save model
        try
            model_save_path = fullfile('results','trainedClassifier_latest.mat');
            trainedModelPackage = struct();
            trainedModelPackage.trainedClassifier = trainedClassifier;
            if isfield(trainedClassifier, 'RequiredVariables')
                trainedModelPackage.requiredVariables = trainedClassifier.RequiredVariables;
            else
                trainedModelPackage.requiredVariables = {};
            end
            trainedModelPackage.validationAccuracy = validationAccuracy;
            if exist('accuracy','var')
                trainedModelPackage.testAccuracy = accuracy/100;
            else
                trainedModelPackage.testAccuracy = NaN;
            end
            trainedModelPackage.savedAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
            try
                trainedModelPackage.compactEnsemble = compact(trainedClassifier.ClassificationEnsemble);
            catch
            end
            save(model_save_path, 'trainedModelPackage');
            ts_stamp = char(datetime("now", "Format", "yyyyMMdd_HHmmss"));
            ts_name = sprintf('trainedClassifier_%s.mat', ts_stamp);
            save(fullfile('results', ts_name), 'trainedModelPackage');
            fprintf('✓ Saved trained classifier to: %s\n', model_save_path);
        catch ME
            fprintf('⚠ Failed to save trained classifier: %s\n', ME.message);
        end
    else
        fprintf('Training feature table empty, skip training.\n');
    end

    fprintf('\nDone (SHHS1 + rpoints pipeline).\n');
catch ME_top
    fprintf('Top-level exception: %s\n', ME_top.message);
end
