function [featureTable, featureNames] = extractHeartbeatFeatures(beatInfo, fs)
% File: extractHeartbeatFeatures.m
% Type: Function
% Usage:
%   [featureTable, featureNames] = extractHeartbeatFeatures(beatInfo, fs)
%
% Description:
%   Extract RR/R/Q/S, QRS duration and area features from beatInfo, build a feature table,
%   normalize numeric columns column-wise, and keep the BeatType label.
%
% Inputs:
%   - beatInfo (struct[]), fs (double)
%
% Outputs:
%   - featureTable (table), featureNames (cellstr)
%
% Dependencies:
%   None (built-ins only)
%
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

featureNames = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', ...
                'QRS_Duration', 'QRS_Area', 'BeatType'}; % last column is label

numBeats = length(beatInfo);
if numBeats == 0
    fprintf('Warning: beatInfo is empty. Cannot extract features.\n');
    % Return empty table with same columns
    featureTable = array2table(zeros(0, length(featureNames)), 'VariableNames', featureNames);
    % Ensure BeatType categorical if possible
    if any(strcmp(featureNames, 'BeatType'))
        featureTable.BeatType = categorical(featureTable.BeatType);
    end
    return;
end

% Initialize feature matrix and labels
% Numeric feature matrix (excluding BeatType)
numericalFeatureNames = setdiff(featureNames, {'BeatType'});
features = NaN(numBeats, length(numericalFeatureNames));

% Handle labels separately
beatTypeLabels = cell(numBeats, 1);

% Compute global R indices for RR intervals
allRIndicesInECG = zeros(numBeats, 1);
for k = 1:numBeats
    % beatInfo(k).rIndex is index within segment
    % beatInfo(k).segmentStartIndex is segment start index in original ECG
    allRIndicesInECG(k) = beatInfo(k).segmentStartIndex + beatInfo(k).rIndex - 1;
end

for i = 1:numBeats
    currentBeat = beatInfo(i);
    segment = currentBeat.segment;
    rIdx = currentBeat.rIndex;
    qIdx = currentBeat.qIndex;
    sIdx = currentBeat.sIndex;
    
    % 1. RR intervals (RR_Prev, RR_Post) in seconds
    if i > 1
        rr_prev = (allRIndicesInECG(i) - allRIndicesInECG(i-1)) / fs;
        features(i, strcmp(numericalFeatureNames, 'RR_Prev')) = rr_prev;
    else
        features(i, strcmp(numericalFeatureNames, 'RR_Prev')) = NaN; % first beat has no previous RR
    end
    
    if i < numBeats
        rr_post = (allRIndicesInECG(i+1) - allRIndicesInECG(i)) / fs;
        features(i, strcmp(numericalFeatureNames, 'RR_Post')) = rr_post;
    else
        features(i, strcmp(numericalFeatureNames, 'RR_Post')) = NaN; % last beat has no next RR
    end
    
    % 2. Amplitude features (R, Q, S) in mV (assuming original ECG is mV)
    if ~isnan(rIdx) && rIdx > 0 && rIdx <= length(segment)
        features(i, strcmp(numericalFeatureNames, 'R_Amplitude')) = segment(rIdx);
    end
    if ~isnan(qIdx) && qIdx > 0 && qIdx <= length(segment)
        features(i, strcmp(numericalFeatureNames, 'Q_Amplitude')) = segment(qIdx);
    end
    if ~isnan(sIdx) && sIdx > 0 && sIdx <= length(segment)
        features(i, strcmp(numericalFeatureNames, 'S_Amplitude')) = segment(sIdx);
    end
    
    % 3. Interval features: QRS duration prefers qrsOnIndex/qrsOffIndex (segment),
    %    otherwise fallback to Q/S/R logic
    qrsOnIdx = NaN; qrsOffIdx = NaN;
    if isfield(currentBeat,'qrsOnIndex'); qrsOnIdx = currentBeat.qrsOnIndex; end
    if isfield(currentBeat,'qrsOffIndex'); qrsOffIdx = currentBeat.qrsOffIndex; end
    if ~isnan(qrsOnIdx) && ~isnan(qrsOffIdx) && qrsOffIdx > qrsOnIdx
        features(i, strcmp(numericalFeatureNames, 'QRS_Duration')) = (qrsOffIdx - qrsOnIdx) / fs;
    elseif ~isnan(qIdx) && ~isnan(sIdx) && sIdx > qIdx
        features(i, strcmp(numericalFeatureNames, 'QRS_Duration')) = (sIdx - qIdx) / fs;
    elseif ~isnan(rIdx) && ~isnan(sIdx) && sIdx > rIdx
        features(i, strcmp(numericalFeatureNames, 'QRS_Duration')) = (sIdx - rIdx) / fs;
    elseif ~isnan(qIdx) && ~isnan(rIdx) && rIdx > qIdx
        features(i, strcmp(numericalFeatureNames, 'QRS_Duration')) = (rIdx - qIdx) / fs;
    else
        features(i, strcmp(numericalFeatureNames, 'QRS_Duration')) = 0.08;
    end
    
    % 4. Area features: QRS area via trapezoidal integration, prefer qrsOn/qrsOff
    if ~isnan(qrsOnIdx) && ~isnan(qrsOffIdx) && qrsOffIdx > qrsOnIdx && qrsOffIdx <= length(segment)
        qrs_segment = segment(qrsOnIdx:qrsOffIdx);
        features(i, strcmp(numericalFeatureNames, 'QRS_Area')) = trapz(abs(qrs_segment)) / fs;
    elseif ~isnan(qIdx) && ~isnan(sIdx) && sIdx > qIdx && sIdx <= length(segment)
        qrs_segment = segment(qIdx:sIdx);
        features(i, strcmp(numericalFeatureNames, 'QRS_Area')) = trapz(abs(qrs_segment)) / fs;
    elseif ~isnan(rIdx)
        r_area_start = max(1, rIdx - round(0.04*fs));
        r_area_end = min(length(segment), rIdx + round(0.04*fs));
        if r_area_end > r_area_start
            qrs_segment = segment(r_area_start:r_area_end);
            features(i, strcmp(numericalFeatureNames, 'QRS_Area')) = trapz(abs(qrs_segment)) / fs;
        end
    end
    
    % No T-related area
    
    % 5. Beat type label
    if isa(currentBeat.beatType, 'char') || isa(currentBeat.beatType, 'string')
        beatTypeLabels{i} = char(currentBeat.beatType); % ensure char
    else
        beatTypeLabels{i} = 'Other'; % default
    end
end

% Create numeric feature table first
numericalFeatureTable = array2table(features, 'VariableNames', numericalFeatureNames);

% Add BeatType column
numericalFeatureTable.BeatType = beatTypeLabels;

% Reorder columns to ensure BeatType is last
featureTable = numericalFeatureTable(:, featureNames);

% Normalize features (excluding BeatType)
fprintf('Start feature normalization...\n');
featureColumns = setdiff(featureNames, {'BeatType'}); % exclude BeatType

for i = 1:length(featureColumns)
    colName = featureColumns{i};
    colData = featureTable{:, colName};
    
    % Remove NaNs for statistics
    validData = colData(~isnan(colData));
    
    if ~isempty(validData) && length(validData) > 1
        % Mean and std
        meanVal = mean(validData);
        stdVal = std(validData);
        
        % Z-score normalization (if std > 0)
        if stdVal > 1e-10 % avoid division by zero
            normalizedData = (colData - meanVal) / stdVal;
            featureTable{:, colName} = normalizedData;
            % fprintf('Feature %s: mean=%.4f, std=%.4f, normalized\n', colName, meanVal, stdVal);
        else
            % If std==0, do min-max normalization
            minVal = min(validData);
            maxVal = max(validData);
            if maxVal > minVal
                normalizedData = (colData - minVal) / (maxVal - minVal);
                featureTable{:, colName} = normalizedData;
                % fprintf('Feature %s: std==0, min-max normalization [%.4f, %.4f]\n', colName, minVal, maxVal);
            else
                fprintf('Feature %s: all values equal, skip normalization\n', colName);
            end
        end
    else
        fprintf('Feature %s: insufficient valid data, skip normalization\n', colName);
    end
end

% Ensure BeatType is string type
if any(strcmp(featureNames, 'BeatType')) && ~isempty(featureTable)
    % Convert BeatType to string array if it is cell
    if iscell(featureTable.BeatType)
        featureTable.BeatType = string(featureTable.BeatType);
    end
    fprintf('BeatType column ensured as string\n');
end

% Feature extraction report
fprintf('Feature extraction finished: %d beats, %d features (normalized).\n', numBeats, length(featureNames)-1);

% Show per-beat-type counts
if ~isempty(featureTable)
    beatTypes = featureTable.BeatType;
    uniqueTypes = unique(beatTypes);
    fprintf('Beat type distribution:\n');
    for i = 1:length(uniqueTypes)
        count = sum(strcmp(beatTypes, uniqueTypes{i}));
        fprintf('  %s: %d beats (%.1f%%)\n', uniqueTypes{i}, count, (count/numBeats)*100);
    end
end

end

