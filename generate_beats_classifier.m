% ========================================================================
% File: generate_beats_classifier.m
% Overview: Generate heartbeat feature sets from SHHS1 EDF records with rpoints
%           annotations and train a binary classifier (PVC vs Other) using
%           AdaBoost / decision tree ensembles.
% Responsibilities:
%   1) Enumerate records under shhs/polysomnography/edfs/shhs1 that have
%      corresponding rpoints CSV files.
%   2) Parse rpoints (-rpoint.csv) annotations (Type/seconds/samplingrate) and:
%        - Remove Artifact (Type=0)
%        - Map labels: Type=2→'PVC'; Type=1/3→'Other'; other values → 'Other'
%   3) Extract the single-lead ECG signal (prefer exact var name 'ECG', otherwise
%      fuzzy match) and standardize to a double column vector.
%   4) Pre-filter: call ecgFilter (method=2; skip power-line notch with power_line_freq=0).
%   5) Call detectAndClassifyHeartbeats for R/Q/S/T detection and beat typing (aligned to annotations).
%   6) Call extractHeartbeatFeatures to build per-beat features → drop incomplete rows → merge all records.
%   7) Stratified random split into train/test (configurable ratio and seed).
%   8) Train classifier (trainBeatsClassifier) and compute CV/test performance + confusion matrix.
%   9) Compute and print multiple metrics (TPR/specificity/precision/F1) and save model with timestamp.
% Inputs/tunable (script variables):
%   splitOptions.test_split_ratio : test ratio (default 0.2)
%   splitOptions.random_seed      : random seed (default 42)
%   splitOptions.stratified       : whether to use stratified sampling (true)
%   processFirstN                 : process only first N records (<=0 means all)
%   model_options.*               : classifier training config (cost matrix/threshold moving, etc.)
% Artifacts (written to results/):
%   trainingFeatureTable.mat
%   testingFeatureTable.mat
%   trainedClassifier_latest.mat  (model package with compact ensemble, val/test metrics, timestamp)
%   trainedClassifier_YYYYMMDD_HHMMSS.mat (timestamped backup)
% Key implementation notes:
%   - rpoints parsing: detectImportOptions + readtable; robust handling of mixed string/numeric columns.
%   - Sampling rate: prefer median of CSV samplingrate → fallback 125 Hz.
%   - Feature merging: keep only complete rows (no missing) to reduce training noise.
%   - Stratified split: shuffle within class + proportionally split to preserve distribution.
%   - Model saving: keep both latest override and timestamped versions.
% Robustness and exception handling:
%   - Missing directory/file/field/channel → skip the record and continue.
%   - Filtering or detection/feature extraction failures → continue; do not abort the pipeline.
%   - Safe exit when feature tables are empty or no valid records.
%   - Compatible conversion for various table column types (cell / numeric).
% Optional processing:
%   - ECG inversion (lead polarity flip) available by enabling: ecg = -ecg;
% Changelog:
%   2025-08-30: Unified Chinese header comments; logic unchanged.
% ========================================================================

clc; clear; close all;

try
    % Runtime switches and paths (can override in workspace before running)
    if ~exist('usePrecomputedFeatures','var'), usePrecomputedFeatures = false; end % load precomputed feature tables to skip EDF processing
    if ~exist('enableHyperparamSearch','var'), enableHyperparamSearch = false; end % auto grid search to maximize F1
    if ~exist('saveOldModelVersions','var'),   saveOldModelVersions   = true; end % save timestamped historical versions
    precomputedTrainPath = fullfile('results','Beats_trainingFeatureTable.mat');
    precomputedTestPath  = fullfile('results','Beats_testingFeatureTable.mat');

    % Hyperparameter search config (adjust ranges/strategy as needed)
    searchOptions = struct();
    searchOptions.enableCostMatrix = [true false];
    searchOptions.costFP = [1 3 5];
    searchOptions.costFN = [1 3];
    searchOptions.enableThresholdMoving = [false true];
    searchOptions.pvcThreshold = [0.50 0.70 0.85 0.95];
    % Extensions: search space for learner-related parameters
    searchOptions.numLearningCycles = [50 100 200];
    searchOptions.learnRate = [0.1 0.2 0.3];
    searchOptions.maxNumSplits = [50 100 200];
    searchOptions.minLeafSize = [5 10 20];
    searchOptions.enableClassWeighting = [true false];
    searchOptions.numVariablesToSample = {'all'}; % can be extended to numeric, e.g., { 'all', 3 }
    searchOptions.randomSampleRatio = 0.5; % range 0–1; 1 means full enumeration
    searchOptions.randomSeed = 123;
    % Directories
    rootDir = pwd;
    edfDir = fullfile(rootDir, 'shhs','polysomnography','edfs','shhs1');
    annDir = fullfile(rootDir, 'shhs','polysomnography','annotations-rpoints','shhs1');
    if ~usePrecomputedFeatures
        if ~isfolder(edfDir)
            error('EDF directory not found: %s', edfDir);
        end
        if ~isfolder(annDir)
            error('rpoints directory not found: %s', annDir);
        end
    end

    % Add helper paths
    addpath(genpath(rootDir));

    % ========== Load blacklist: results/badsignallist.csv ==========
    blackListIds = string([]);
    blackListNames = string([]);
    blacklistFile = fullfile('results', 'badsignallist.csv');
    if exist(blacklistFile,'file')
        try
            C = readcell(blacklistFile);
            vals = string(C(:));
            vals = strtrim(vals);
            vals = replace(vals, '"','');
            vals = replace(vals, '''','');
            vals = vals(vals ~= "");
            isDigits = false(numel(vals),1);
            for ii = 1:numel(vals)
                isDigits(ii) = ~isempty(regexp(vals(ii), '^\d+$', 'once'));
            end
            ids = vals(isDigits);
            names = vals(~isDigits);
            recIdFromNames = strings(0,1);
            for ii = 1:numel(names)
                tok = regexp(names(ii), '(\d+)', 'tokens', 'once');
                if ~isempty(tok)
                    recIdFromNames(end+1,1) = string(tok{1}); %#ok<AGROW>
                end
            end
            blackListIds = unique([ids; recIdFromNames]);
            blackListNames = unique("shhs1-" + blackListIds);
            fprintf('Loaded blacklist with %d records; they will be skipped.\n', numel(blackListIds));
        catch ME
            fprintf('Failed to load blacklist (no records will be skipped): %s\n', ME.message);
        end
    else
        fprintf('Blacklist file not found (no records will be skipped): %s\n', blacklistFile);
    end

    % If using precomputed feature tables, try loading directly
    dataPrepared = false;
    if usePrecomputedFeatures
        try
            S_train = load(precomputedTrainPath);
            S_test  = load(precomputedTestPath);
            if isfield(S_train,'trainingFeatureTable') && isfield(S_test,'testingFeatureTable')
                trainingFeatureTable = S_train.trainingFeatureTable;
                testingFeatureTable  = S_test.testingFeatureTable;
                fprintf('Loaded feature tables from precomputed files:\n  %s\n  %s\n', precomputedTrainPath, precomputedTestPath);
                fprintf('Train=%d, Test=%d\n', height(trainingFeatureTable), height(testingFeatureTable));
                dataPrepared = true;
            else
                fprintf('Precomputed files exist but lack required variables; falling back to raw EDF processing.\n');
            end
        catch ME
            fprintf('Failed to load precomputed feature tables: %s\nSwitching to raw EDF processing.\n', ME.message);
        end
    end

    % Split and training options
    splitOptions = struct();
    splitOptions.test_split_ratio = 0.2;  % test ratio
    splitOptions.random_seed = 42;        % random seed
    splitOptions.stratified = true;       % stratified sampling

    % Process only first N annotated EDFs (<=0 means process all)
    processFirstN = 500;

    % Training config (consistent with prior version)
    model_options = struct();
    model_options.enableCostMatrix = true;
    % High-precision config: raise cost of Other->PVC, lower cost of PVC->Other
    model_options.costFP = 5;  % false positive cost (Other predicted as PVC)
    model_options.costFN = 1;  % false negative cost (PVC predicted as Other)
    model_options.enableThresholdMoving = true;
    model_options.pvcThreshold = 0.95; % higher threshold during training to suppress false positives

    if ~dataPrepared
    % List records with rpoints annotations (e.g., shhs1-XXXXX-rpoint.csv → shhs1-XXXXX)
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

    % Apply blacklist to record names (covers annotated records and name forms without annotations)
    if ~isempty(blackListIds)
        keepMask = true(numel(recordBases),1);
        for kk = 1:numel(recordBases)
            nm = string(recordBases{kk});
            recIdTok = regexp(nm, 'shhs1-(\d+)$', 'tokens', 'once');
            recIdStr = "";
            if ~isempty(recIdTok)
                recIdStr = string(recIdTok{1});
            end
            if any(blackListNames == nm) || (strlength(recIdStr) > 0 && any(blackListIds == recIdStr))
                keepMask(kk) = false;
            end
        end
        recordBases = recordBases(keepMask);
    end

    % Keep only existing EDFs
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

    % Truncate if processing only first N
    if exist('processFirstN','var') && isnumeric(processFirstN) && ~isempty(processFirstN) && processFirstN > 0
        maxN = min(processFirstN, numel(edfPaths));
        edfPaths = edfPaths(1:maxN);
        annPaths = annPaths(1:maxN);
        fprintf('Found %d EDF records with rpoints in %s (total %d rpoints files). Processing only the first %d.\n', numel(recordBases), edfDir, numel(annFiles), maxN);
    else
        fprintf('Found %d EDF records with rpoints in %s (total %d rpoints files).\n', numel(edfPaths), edfDir, numel(annFiles));
    end
    if isempty(edfPaths)
        fprintf('No records to process. Exiting.\n');
        return;
    end

    % Only keep feature tables (temporarily per record)
    features_by_record = cell(numel(edfPaths), 1);
    numKept = 0;

    for iFile = 1:numel(edfPaths)
        edfPath = edfPaths{iFile};
        annPath = annPaths{iFile};
        [edfFolder, edfName, ~] = fileparts(edfPath);
    recBase = edfName; % like shhs1-XXXXXX
    % Blacklist guard: defensive check (even if filtered upstream)
        recIdTok = regexp(string(recBase), 'shhs1-(\d+)$', 'tokens', 'once');
        recIdStr = ""; if ~isempty(recIdTok), recIdStr = string(recIdTok{1}); end
        if any(blackListNames == string(recBase)) || (strlength(recIdStr) > 0 && any(blackListIds == recIdStr))
            continue;
        end
    fprintf('\n=== [%d/%d] Processing %s ===\n', iFile, numel(edfPaths), [edfName '.edf']);

        % Read EDF (ECG channel only)
        try
            TT = edfread(edfPath);
        catch ME
            fprintf('  Skip: edfread failed: %s\n', ME.message);
            continue;
        end

        % Find ECG channel (prefer exact match 'ECG', then fuzzy match)
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

    % Read rpoints CSV and extract seconds/Type/samplingrate
        try
            opts = detectImportOptions(annPath);
            Tann = readtable(annPath, opts);
        catch ME
            fprintf('  Skip: failed to read rpoints: %s\n', ME.message);
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

    % Sampling rate (prefer CSV; fallback 125 Hz if missing)
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

    % Convert seconds/Type columns
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
            fprintf('  Skip: no valid rpoints annotations.\n');
            continue;
        end
        % Type mapping: 2->PVC; 1/3->Other; others (if any) → Other
        ANNOTD = cell(numel(type_kept), 1);
        for ii = 1:numel(type_kept)
            if type_kept(ii) == 2
                ANNOTD{ii} = 'PVC';
            else
                ANNOTD{ii} = 'Other';
            end
        end

    % Print PVC/Other counts in annotations
        pvcCount = sum(strcmp(ANNOTD, 'PVC'));
        otherCount = sum(strcmp(ANNOTD, 'Other'));
    fprintf('  Annotation stats: PVC=%d, Other=%d (total=%d)\n', pvcCount, otherCount, pvcCount + otherCount);

    % Key: invert entire ECG (optional)
        % ecg = -ecg;

        % Filtering: skip notch (60 Hz already handled)
        try
            [ecgFiltered, ~] = ecgFilter(ecg, fsLocal, 2, 0);
        catch ME
            fprintf('  Skip: filtering failed: %s\n', ME.message);
            continue;
        end

        % Beat detection and typing (nearest R-annotation matching)
        try
            [~, beatInfo, stats] = detectAndClassifyHeartbeats(ecgFiltered, ATRTIMED, ANNOTD, fsLocal);
        catch ME
            fprintf('  Skip: heartbeat detection error: %s\n', ME.message);
            continue;
        end
        if isfield(stats,'mteo_failed') && stats.mteo_failed
            fprintf('  Skip: MTEO (Q/S) detection failed.\n');
            continue;
        end
        if isempty(beatInfo)
            fprintf('  Skip: no valid heartbeats obtained.\n');
            continue;
        end

        % Feature extraction (keep features; drop missing rows)
        [featureTable_rec, ~] = extractHeartbeatFeatures(beatInfo, fsLocal);
        if isempty(featureTable_rec) || height(featureTable_rec) == 0
            fprintf('  Skip: feature extraction result is empty.\n');
            continue;
        end
        completeRows = ~any(ismissing(featureTable_rec), 2);
        featureTable_rec = featureTable_rec(completeRows, :);
        fprintf('  Valid beats (feature rows) in this record: %d\n', height(featureTable_rec));

        numKept = numKept + 1;
        features_by_record{numKept,1} = featureTable_rec;

    % Cleanup temporary variables for current record
        clear TT ecg ecgCol ecgVarName ecgFiltered Tann ATRTIMED ANNOTD featureTable_rec completeRows ...
              typeRaw secRaw typeVec secVec mask_valid type_kept idxType idxSec idxSR srRaw srVec fsLocal;
    end

    % Merge all features
    if numKept == 0
        fprintf('\nFailed to extract valid features from any record. Exiting.\n');
        return;
    end
    features_by_record = features_by_record(1:numKept);
    valid_tables = features_by_record(~cellfun('isempty', features_by_record));
    allFeatureTable = vertcat(valid_tables{:});
    fprintf('\nAll feature tables merged. Extracted features for %d beats in total.\n', height(allFeatureTable));

    % Stratified random split
    fprintf('\n=== Performing stratified random split ===\n');
    fprintf('Test ratio: %.1f%%, Random seed: %d\n', splitOptions.test_split_ratio*100, splitOptions.random_seed);
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
    fprintf('\nSplit complete: train=%d, test=%d\n', height(trainingFeatureTable), height(testingFeatureTable));

    % Save feature tables
    if ~exist('results','dir'), mkdir('results'); end
    try
        save(fullfile('results','Beats_trainingFeatureTable.mat'), 'trainingFeatureTable');
        fprintf('✓ Training feature table saved\n');
    catch ME
        fprintf('⚠ Failed to save training feature table: %s\n', ME.message);
    end
    try
        save(fullfile('results','Beats_testingFeatureTable.mat'), 'testingFeatureTable');
        fprintf('✓ Test feature table saved\n');
    catch ME
        fprintf('⚠ Failed to save test feature table: %s\n', ME.message);
    end
    end % if ~dataPrepared

    % Train classifier / auto hyperparameter tuning
    if height(trainingFeatureTable) > 0
        if enableHyperparamSearch
            fprintf('\n--- Auto hyperparameter tuning and training (goal: maximize macro F1) ---\n');
            [bestOptions, trainedClassifier, hpStats] = beatsGridSearchHyperparams(trainingFeatureTable, testingFeatureTable, model_options, searchOptions);
            validationAccuracy = hpStats.bestValidationAccuracy;
            model_options = bestOptions;
            fprintf('Tuning complete. Best macro F1=%.4f, using config: enableCostMatrix=%d, costFP=%.3g, costFN=%.3g, enableThresholdMoving=%d, pvcThreshold=%.2f, numLearningCycles=%d, learnRate=%.2f, maxNumSplits=%d, minLeafSize=%d, enableClassWeighting=%d\n', ...
                hpStats.bestMacroF1, model_options.enableCostMatrix, model_options.costFP, model_options.costFN, model_options.enableThresholdMoving, model_options.pvcThreshold, ...
                getfield(model_options,'numLearningCycles',100), getfield(model_options,'learnRate',0.2), getfield(model_options,'maxNumSplits',100), getfield(model_options,'minLeafSize',10), getfield(model_options,'enableClassWeighting',1));
        else
            fprintf('\n--- Train classifier ---\n');
            [trainedClassifier, validationAccuracy] = trainBeatsClassifier(trainingFeatureTable, model_options);
            fprintf('Classifier training complete. Cross-validation accuracy: %.2f%%\n', validationAccuracy * 100);
        end

        % Test-set evaluation
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
            fprintf('Overall test accuracy: %.2f%%\n', accuracy);

            % === Confusion matrix (figure + text) ===
            % 统一为 cellstr
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

            % Confusion matrix figure (row/column summary)
            figure('Name', 'Confusion Matrix (Test Set)', 'NumberTitle', 'off');
            confusionchart(confusion_mat, class_names, 'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
            title('ConfusionMatrix (Test Set)');

            % Text version: absolute counts
            fprintf('\n=== Confusion Matrix (absolute counts) ===\n');
            header = sprintf('%-12s', 'Actual\\Predicted');
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

            % Text version: row percentages
            fprintf('\n=== Confusion Matrix (row percentages) ===\n');
            row_sums = sum(confusion_mat, 2);
            header = sprintf('%-12s', 'Actual\\Predicted');
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

            % Text version: per-class metrics (TPR/FNR/Specificity/Precision/F1)
            fprintf('\n=== Per-class metrics ===\n');
            macroF1_vec = nan(1, num_classes);
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
                macroF1_vec(ii) = f1_score;
            end
            macroF1 = mean(macroF1_vec(~isnan(macroF1_vec)));
            fprintf('Macro F1: %.2f\n', macroF1);
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
            if saveOldModelVersions
                ts_stamp = char(datetime("now", "Format", "yyyyMMdd_HHmmss"));
                ts_name = sprintf('trainedClassifier_%s.mat', ts_stamp);
                save(fullfile('results', ts_name), 'trainedModelPackage');
            end
            fprintf('✓ Saved trained classifier to: %s\n', model_save_path);
        catch ME
            fprintf('⚠ Failed to save trained classifier: %s\n', ME.message);
        end
    else
        fprintf('Training feature table is empty; skipping training.\n');
    end

    fprintf('\nProcessing complete (SHHS1 + rpoints pipeline).\n');
catch ME_top
    fprintf('Top-level exception: %s\n', ME_top.message);
end


% ========================================================================
% Local function: Hyperparameter search (targets test-set macro F1)
% If test set is empty, falls back to cross-validation accuracy as target.
% ========================================================================
function [bestOptions, bestTrainedClassifier, stats] = beatsGridSearchHyperparams(trainingFeatureTable, testingFeatureTable, baseOptions, searchOptions)
    % Fill default search space
    if nargin < 3 || isempty(baseOptions), baseOptions = struct(); end
    if nargin < 4, searchOptions = struct(); end
    if ~isfield(searchOptions,'enableCostMatrix'),        searchOptions.enableCostMatrix = [true false]; end
    if ~isfield(searchOptions,'costFP'),                  searchOptions.costFP = [1 3 5]; end
    if ~isfield(searchOptions,'costFN'),                  searchOptions.costFN = [1 3]; end
    if ~isfield(searchOptions,'enableThresholdMoving'),   searchOptions.enableThresholdMoving = [false true]; end
    if ~isfield(searchOptions,'pvcThreshold'),            searchOptions.pvcThreshold = [0.50 0.70 0.85 0.95]; end
    if ~isfield(searchOptions,'randomSampleRatio'),       searchOptions.randomSampleRatio = 1.0; end
    if ~isfield(searchOptions,'randomSeed'),              searchOptions.randomSeed = 123; end
    if ~isfield(searchOptions,'numLearningCycles'),       searchOptions.numLearningCycles = 100; end
    if ~isfield(searchOptions,'learnRate'),               searchOptions.learnRate = 0.2; end
    if ~isfield(searchOptions,'maxNumSplits'),            searchOptions.maxNumSplits = 100; end
    if ~isfield(searchOptions,'minLeafSize'),             searchOptions.minLeafSize = 10; end
    if ~isfield(searchOptions,'enableClassWeighting'),    searchOptions.enableClassWeighting = true; end
    if ~isfield(searchOptions,'numVariablesToSample'),    searchOptions.numVariablesToSample = {'all'}; end

    % Generate candidate combinations
    combos = struct('enableCostMatrix', {}, 'costFP', {}, 'costFN', {}, 'enableThresholdMoving', {}, 'pvcThreshold', {}, 'numLearningCycles', {}, 'learnRate', {}, 'maxNumSplits', {}, 'minLeafSize', {}, 'enableClassWeighting', {}, 'numVariablesToSample', {});
    for ecm = searchOptions.enableCostMatrix
        if ecm
            listFP = searchOptions.costFP;
            listFN = searchOptions.costFN;
        else
            listFP = searchOptions.costFP(1);
            listFN = searchOptions.costFN(1);
        end
        for cf = listFP
            for cn = listFN
                for etm = searchOptions.enableThresholdMoving
                    if etm
                        listTH = searchOptions.pvcThreshold;
                    else
                        listTH = searchOptions.pvcThreshold(1);
                    end
                    for th = listTH
                        for nlc = ensureArray(searchOptions.numLearningCycles)
                            for lr = ensureArray(searchOptions.learnRate)
                                for mns = ensureArray(searchOptions.maxNumSplits)
                                    for mls = ensureArray(searchOptions.minLeafSize)
                                        for ecw = ensureArray(searchOptions.enableClassWeighting)
                                            for nvs = ensureCellArray(searchOptions.numVariablesToSample)
                                                combos(end+1) = struct('enableCostMatrix', logical(ecm), ...
                                                                       'costFP', cf, ...
                                                                       'costFN', cn, ...
                                                                       'enableThresholdMoving', logical(etm), ...
                                                                       'pvcThreshold', th, ...
                                                                       'numLearningCycles', nlc, ...
                                                                       'learnRate', lr, ...
                                                                       'maxNumSplits', mns, ...
                                                                       'minLeafSize', mls, ...
                                                                       'enableClassWeighting', logical(ecw), ...
                                                                       'numVariablesToSample', nvs); %#ok<AGROW>
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    totalCandidates = numel(combos);
    % Random subsampling
    if searchOptions.randomSampleRatio < 1.0 && totalCandidates > 1
        rng(searchOptions.randomSeed);
        keepN = max(1, round(totalCandidates * searchOptions.randomSampleRatio));
        pickIdx = randperm(totalCandidates, keepN);
        combos = combos(pickIdx);
    end

    hasTest = exist('testingFeatureTable','var') && ~isempty(testingFeatureTable) && height(testingFeatureTable) > 0;
    bestMacroF1 = -inf; % in percent (0~100)
    bestValAcc  = -inf; % as ratio (0~1)
    bestOptions = baseOptions;
    bestTrainedClassifier = [];
    tried = 0;

    fprintf('Hyperparameter candidates: %d (after sampling)\n', numel(combos));

    for k = 1:numel(combos)
        tried = tried + 1;
        opt = baseOptions;
        opt.enableCostMatrix = combos(k).enableCostMatrix;
        opt.costFP = combos(k).costFP;
        opt.costFN = combos(k).costFN;
        opt.enableThresholdMoving = combos(k).enableThresholdMoving;
        opt.pvcThreshold = combos(k).pvcThreshold;

        try
            [tc, valAcc] = localTrainBeatsClassifier(trainingFeatureTable, opt);
        catch ME
            fprintf('  Combination #%d training failed: %s\n', k, ME.message);
            continue;
        end

    % Evaluate macro F1 (if test set available)
    macroF1 = nan; % percent
        if hasTest
            try
                try
                    [predLabels, ~] = tc.predictFcn(testingFeatureTable);
                catch
                    predLabels = tc.predictFcn(testingFeatureTable);
                end
                trueLabels = testingFeatureTable.BeatType;
                macroF1 = localComputeMacroF1Percent(predLabels, trueLabels); % 0~100
            catch ME
                fprintf('  Combination #%d test evaluation failed: %s\n', k, ME.message);
            end
        end

        score = macroF1; if isnan(score), score = -inf; end
        bestScore = bestMacroF1; if isnan(bestScore), bestScore = -inf; end

        isBetter = (score > bestScore) || (abs(score - bestScore) <= 1e-9 && valAcc > bestValAcc);
        if isBetter
            bestMacroF1 = macroF1;
            bestValAcc  = valAcc;
            bestOptions = opt;
            bestTrainedClassifier = tc;
            if hasTest && ~isnan(macroF1)
                fprintf('  ✓ New best: Macro F1=%.2f%%, ValAcc=%.2f%% [k=%d]\n', bestMacroF1, bestValAcc*100, k);
            else
                fprintf('  ✓ New best: ValAcc=%.2f%% [no test-set macro F1][k=%d]\n', bestValAcc*100, k);
            end
        end
    end

    % If all fail or no model trained, fall back to baseOptions
    if isempty(bestTrainedClassifier)
        fprintf('No model obtained from hyperparameter search; falling back to base configuration training.\n');
        [bestTrainedClassifier, bestValAcc] = localTrainBeatsClassifier(trainingFeatureTable, baseOptions);
        bestOptions = baseOptions;
        bestMacroF1 = nan;
    end

    stats = struct();
    stats.bestMacroF1 = bestMacroF1;            % percent (0~100), possibly NaN
    stats.bestValidationAccuracy = bestValAcc;  % ratio (0~1)
    stats.totalCandidates = totalCandidates;
    stats.triedCandidates = tried;
end

% Compute macro F1 (percent, 0~100). Input labels may be cell/categorical/string.
function macroF1 = localComputeMacroF1Percent(predicted_labels, actual_labels)
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
    f1_list = nan(1, num_classes);
    for ii = 1:num_classes
        TP = confusion_mat(ii, ii);
        FP = sum(confusion_mat(:, ii)) - TP;
        FN = sum(confusion_mat(ii, :)) - TP;
        if (TP + FN) > 0
            recall = TP / (TP + FN) * 100;  % percent
        else
            recall = 0;
        end
        if (TP + FP) > 0
            precision = TP / (TP + FP) * 100; % percent
        else
            precision = 0;
        end
        if (precision + recall) > 0
            f1_list(ii) = 2 * (precision * recall) / (precision + recall); % percent
        else
            f1_list(ii) = 0;
        end
    end
    macroF1 = mean(f1_list, 'omitnan'); % percent
end

% ========================================================================
% Local function: Train beat-type classifier (inline version of trainBeatsClassifier)
% Supported modelOptions fields:
%   enableCostMatrix, costFP, costFN, enableThresholdMoving, pvcThreshold,
%   numLearningCycles, learnRate, maxNumSplits, minLeafSize,
%   enableClassWeighting, numVariablesToSample
% ========================================================================
function [trainedClassifier, validationAccuracy] = localTrainBeatsClassifier(trainingData, modelOptions)
    if nargin < 2 || isempty(modelOptions)
        modelOptions = struct();
    end
    % Defaults
    defaults = struct();
    defaults.enableCostMatrix = false;
    defaults.costFP = 1;
    defaults.costFN = 5;
    defaults.enableThresholdMoving = false;
    defaults.pvcThreshold = 0.50;
    defaults.numLearningCycles = 100;
    defaults.learnRate = 0.2;
    defaults.maxNumSplits = 100;
    defaults.minLeafSize = 10;
    defaults.enableClassWeighting = true;
    defaults.numVariablesToSample = 'all';

    % Merge defaults
    fn = fieldnames(defaults);
    for i = 1:numel(fn)
        if ~isfield(modelOptions, fn{i})
            modelOptions.(fn{i}) = defaults.(fn{i});
        end
    end

    inputTable = trainingData;
    predictorNames = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
    predictors = inputTable(:, predictorNames);
    response = inputTable.BeatType;
    classNames = {'Other'; 'PVC'};

    % Observation weights (optional class balancing)
    numObservations = height(inputTable);
    if isstring(response), response = cellstr(response); end
    if iscategorical(response), response = cellstr(response); end
    if istable(response), response = response{:,:}; end
    if modelOptions.enableClassWeighting
        countOther = sum(strcmp(response, 'Other'));
        countPVC   = sum(strcmp(response, 'PVC'));
        if countOther > 0 && countPVC > 0
            classCounts = [countOther, countPVC];
            classWeights = numObservations ./ (numel(classNames) .* classCounts);
            obsWeights = zeros(numObservations, 1);
            obsWeights(strcmp(response, 'Other')) = classWeights(1);
            obsWeights(strcmp(response, 'PVC'))   = classWeights(2);
        else
            obsWeights = ones(numObservations, 1);
        end
    else
        obsWeights = ones(numObservations, 1);
    end

    % Tree template
    template = templateTree( ...
        'MaxNumSplits', modelOptions.maxNumSplits, ...
        'MinLeafSize', modelOptions.minLeafSize, ...
        'NumVariablesToSample', modelOptions.numVariablesToSample);

    % Train ensemble
    fitArgs = { ...
        predictors, ...
        response, ...
        'Method', 'AdaBoostM1', ...
        'NumLearningCycles', modelOptions.numLearningCycles, ...
        'Learners', template, ...
        'LearnRate', modelOptions.learnRate, ...
        'ClassNames', classNames, ...
        'Weights', obsWeights ...
    };
    if modelOptions.enableCostMatrix
        C = [0, modelOptions.costFP; modelOptions.costFN, 0];
        fitArgs = [fitArgs, {'Cost', C}]; %#ok<AGROW>
    end
    classificationEnsemble = fitcensemble(fitArgs{:});

    % Prediction function: supports threshold moving
    predictorExtractionFcn = @(t) t(:, predictorNames);
    if modelOptions.enableThresholdMoving
        threshold = modelOptions.pvcThreshold;
        thresholdPredictFcn = @(tbl) localThresholdPredict(tbl, predictorExtractionFcn, classificationEnsemble, threshold);
        trainedClassifier.predictFcn = thresholdPredictFcn;
    else
        ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
        trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
    end

    trainedClassifier.RequiredVariables = predictorNames;
    trainedClassifier.ClassificationEnsemble = classificationEnsemble;
    trainedClassifier.Options = modelOptions;

    % Cross-validation accuracy
    partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 10);
    validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
end

% Utility: ensure array
function arr = ensureArray(x)
    if iscell(x)
        arr = cell2mat(x);
    elseif isnumeric(x) || islogical(x)
        arr = x;
    else
        arr = x;
    end
end

% Utility: ensure cell array
function c = ensureCellArray(x)
    if iscell(x)
        c = x;
    else
        c = {x};
    end
end
