function [trainedClassifier, validationAccuracy] = trainBeatsClassifier(trainingData, modelOptions)
% ========================================================================
% File: trainBeatsClassifier.m
% Overview: Train a beat-type classifier (PVC vs Other) from a feature table using
%           an AdaBoostM1 ensemble of decision trees. Supports class-frequency
%           observation weights, an optional cost matrix, and optional threshold
%           moving to trade PVC precision/recall.
% Responsibilities:
%   1) Accept a table (trainingData) with numeric predictors and a BeatType
%      label ('Other'/'PVC').
%   2) Assemble the predictor set (RR / amplitudes / QRS duration / area) and
%      build observation weights for class balancing.
%   3) Optionally apply a cost matrix (false positive costFP, false negative
%      costFN) to adjust the decision boundary.
%   4) Train an AdaBoostM1 ensemble (tree template: MaxNumSplits=100,
%      MinLeafSize=10, LearnRate=0.2).
%   5) If threshold moving is enabled (enableThresholdMoving), use pvcThreshold
%      to re-threshold the PVC probability/score.
%   6) Return a structured output: the trained model, predictFcn, required
%      variable list, and 10-fold cross-validation accuracy.
% Inputs:
%   trainingData (table):
%       Required columns: {'RR_Prev','RR_Post','R_Amplitude','Q_Amplitude',
%                          'S_Amplitude','QRS_Duration','QRS_Area','BeatType'}
%       BeatType label values: 'Other' or 'PVC'
%   modelOptions (struct, optional):
%       .enableCostMatrix (logical)  : whether to enable a 2x2 cost matrix (default false)
%       .costFP (double)             : cost C(1,2) for Other->PVC (default 1)
%       .costFN (double)             : cost C(2,1) for PVC->Other (default 5)
%       .enableThresholdMoving       : whether to apply post-hoc threshold moving (default false)
%       .pvcThreshold (double)       : threshold in [0,1] (default 0.50)
% Outputs:
%   trainedClassifier (struct):
%       .predictFcn(table) -> [labels,scores] (if threshold moving enabled, returns labels + raw scores)
%       .ClassificationEnsemble        : underlying fitcensemble model
%       .RequiredVariables (cellstr)   : required predictor list
%       .Options                       : modelOptions used for training
%       .About / .HowToPredict         : metadata
%   validationAccuracy (double)        : 10-fold cross-validation accuracy (0-1)
% Key implementation notes:
%   - Class-frequency weighting: inverse-frequency per class to mitigate imbalance.
%   - Cost matrix: affects cost-sensitive decisions during boosting.
%   - Threshold moving: compare rawScores(:,PVC) against pvcThreshold to tune recall/precision.
%   - Cross-validation: crossval(K=10) for generalization performance.
%   - Robust label normalization: string/categorical converted to cellstr.
% Edge cases and robustness:
%   - If one class is missing: fall back to equal weights (avoid divide-by-zero/Inf).
%   - Missing modelOptions fields are filled with defaults.
%   - With threshold moving, preserve the original score matrix for later tuning.
% Changelog:
%   2025-08-30: English header normalization; training logic unchanged.
% ========================================================================


% Extract predictors and response
% The following code prepares data in the proper shape for training.
%
% The second argument is optional
if nargin < 2 || isempty(modelOptions)
    modelOptions = struct();
end
if ~isfield(modelOptions, 'enableCostMatrix'),      modelOptions.enableCostMatrix = false; end
if ~isfield(modelOptions, 'costFP'),                modelOptions.costFP = 1; end
if ~isfield(modelOptions, 'costFN'),                modelOptions.costFN = 5; end
if ~isfield(modelOptions, 'enableThresholdMoving'), modelOptions.enableThresholdMoving = false; end
if ~isfield(modelOptions, 'pvcThreshold'),          modelOptions.pvcThreshold = 0.50; end

inputTable = trainingData;
predictorNames = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
predictors = inputTable(:, predictorNames);
response = inputTable.BeatType;
isCategoricalPredictor = [false, false, false, false, false, false, false]; %#ok<NASGU>
classNames = {'Other'; 'PVC'};

% Compute observation weights from class frequencies (without discarding samples)
numObservations = height(inputTable);
numClasses = numel(classNames);
% Normalize response to cellstr for consistent comparisons
if isstring(response)
    response = cellstr(response);
end
if iscategorical(response)
    response = cellstr(response);
end
if istable(response)
    response = response{:,:};
end
countOther = sum(strcmp(response, 'Other'));
countPVC   = sum(strcmp(response, 'PVC'));
if countOther > 0 && countPVC > 0
    classCounts = [countOther, countPVC];
    classWeights = numObservations ./ (numClasses .* classCounts);
    obsWeights = zeros(numObservations, 1);
    obsWeights(strcmp(response, 'Other')) = classWeights(1);
    obsWeights(strcmp(response, 'PVC'))   = classWeights(2);
else
    % If one class is absent in the training set, fall back to equal weights to avoid div-by-zero/Inf
    obsWeights = ones(numObservations, 1);
end

% Train classifier
% The following sets all classifier options and fits the model.
template = templateTree(...
    'MaxNumSplits', 100, ...
    'MinLeafSize', 10, ...
    'NumVariablesToSample', 'all');
% Build optional 2x2 cost matrix, order corresponds to {'Other','PVC'}
fitArgs = {...
    predictors, ...
    response, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles', 100, ...
    'Learners', template, ...
    'LearnRate', 0.2, ...
    'ClassNames', classNames, ...
    'Weights', obsWeights ...
};

if modelOptions.enableCostMatrix
    % Cost matrix C(i,j): cost of predicting j when the true class is i
    % Higher cost for false negatives (PVC->Other): C(2,1)=costFN; false positives (Other->PVC): C(1,2)=costFP
    C = [0, modelOptions.costFP; modelOptions.costFN, 0];
    fitArgs = [fitArgs, {'Cost', C}]; %#ok<AGROW>
end

classificationEnsemble = fitcensemble(fitArgs{:});

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
% Predict function: supports optional threshold moving
if modelOptions.enableThresholdMoving
    % Return both labels and scores
    threshold = modelOptions.pvcThreshold;
    thresholdPredictFcn = @(tbl) localThresholdPredict(tbl, predictorExtractionFcn, classificationEnsemble, threshold);
    trainedClassifier.predictFcn = thresholdPredictFcn;
else
    ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
    trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
end

% Populate fields on the result struct
trainedClassifier.RequiredVariables = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.Options = modelOptions; % Store training-time options and threshold
trainedClassifier.About = 'This struct contains a trained model exported from Classification Learner R2025a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table T, use:\n [yfit,scores] = c.predictFcn(T) \nReplace ''c'' with the variable name of this struct, e.g., ''trainedModel''.\n \nTable T must contain the variables returned by:\n c.RequiredVariables \nThe formats (e.g., matrix/vector, data type) must match the original training data.\nOther variables in T are ignored.\n \nFor details, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response (duplicate block kept as in exported template)
% The following code prepares data in the proper shape for training.
%
inputTable = trainingData;
predictorNames = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
predictors = inputTable(:, predictorNames);
response = inputTable.BeatType;
isCategoricalPredictor = [false, false, false, false, false, false, false]; %#ok<NASGU>
classNames = {'Other'; 'PVC'}; %#ok<NASGU>

% Cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 10);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel); %#ok<NASGU,ASGLU>

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end

function [labels, scores] = localThresholdPredict(tbl, predictorExtractionFcn, model, pvcThreshold)
% Return labels and scores applying a custom threshold. Assumes the PVC score
% resides in the column corresponding to class name 'PVC'. If the score
% columns align with class names, then scores(:, idxPVC) is the PVC score.
    X = predictorExtractionFcn(tbl);
    [rawLabels, rawScores] = predict(model, X);
    % Locate the PVC score column
    classNames = model.ClassNames; % {'Other','PVC'}
    if iscell(classNames)
        pvcIdx = find(strcmp(classNames, 'PVC'), 1);
    else
        pvcIdx = find(classNames == 'PVC', 1);
    end
    if isempty(pvcIdx)
        % Fallback: use the second column by default
        pvcIdx = min(size(rawScores,2), 2);
    end
    pvcScore = rawScores(:, pvcIdx);
    labels = cell(size(rawLabels));
    labels(:) = {'Other'};
    labels(pvcScore >= pvcThreshold) = {'PVC'};
    scores = rawScores;
end
