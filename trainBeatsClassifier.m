function [trainedClassifier, validationAccuracy] = trainBeatsClassifier(trainingData, modelOptions)
% File: trainBeatsClassifier.m
% Type: Function
% Usage:
%   [trainedClassifier, validationAccuracy] = trainBeatsClassifier(trainingData, modelOptions)
%
% Description:
%   Train a binary beat classifier (Other vs PVC) using AdaBoost trees with class-frequency weights;
%   optional cost matrix and threshold moving via modelOptions.
%
% Inputs:
%   - trainingData (table)
%   - modelOptions (struct, optional): enableCostMatrix/costFP/costFN/enableThresholdMoving/pvcThreshold
%
% Outputs:
%   - trainedClassifier (struct) and validationAccuracy (double)
%
% Dependencies:
%   Statistics and Machine Learning Toolbox built-ins
%
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

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
isCategoricalPredictor = [false, false, false, false, false, false, false];
classNames = {'Other'; 'PVC'};

% Compute observation weights from class frequencies (no discard)
numObservations = height(inputTable);
numClasses = numel(classNames);
% Normalize response to cellstr for comparison
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
    % If a class is missing in training set, fallback to equal weights
    obsWeights = ones(numObservations, 1);
end

% Train classifier: set options and train
template = templateTree(...
    'MaxNumSplits', 100, ...
    'MinLeafSize', 10, ...
    'NumVariablesToSample', 'all');
% Optional cost matrix (2x2), order: {'Other','PVC'}
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
    % Cost matrix C(i,j): true i predicted j
    % FN (PVC->Other) cost high: C(2,1)=costFN; FP (Other->PVC) cost: C(1,2)=costFP
    C = [0, modelOptions.costFP; modelOptions.costFN, 0];
    fitArgs = [fitArgs, {'Cost', C}]; %#ok<AGROW>
end

classificationEnsemble = fitcensemble(fitArgs{:});

% Create result struct
predictorExtractionFcn = @(t) t(:, predictorNames);
% Prediction function with threshold moving
if modelOptions.enableThresholdMoving
    threshold = modelOptions.pvcThreshold;
    thresholdPredictFcn = @(tbl) localThresholdPredict(tbl, predictorExtractionFcn, classificationEnsemble, threshold);
    trainedClassifier.predictFcn = thresholdPredictFcn;
else
    ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
    trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
end

% Add fields to result struct
trainedClassifier.RequiredVariables = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.Options = modelOptions;
trainedClassifier.About = 'Model exported from Classification Learner R2025a.';
trainedClassifier.HowToPredict = sprintf('To predict using a new table T:\n [yfit,scores] = c.predictFcn(T)\nReplace ''c'' with the variable name (e.g., ''trainedModel'').\n\nTable T must include variables returned by:\n c.RequiredVariables\nFormats must match the training data; other variables are ignored.\n\nSee <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors/response and prepare data for training
%
inputTable = trainingData;
predictorNames = {'RR_Prev', 'RR_Post', 'R_Amplitude', 'Q_Amplitude', 'S_Amplitude', 'QRS_Duration', 'QRS_Area'};
predictors = inputTable(:, predictorNames);
response = inputTable.BeatType;
isCategoricalPredictor = [false, false, false, false, false, false, false];
classNames = {'Other'; 'PVC'};

% Cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 10);

% Validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end

function [labels, scores] = localThresholdPredict(tbl, predictorExtractionFcn, model, pvcThreshold)
% Return labels/scores with custom threshold; PVC score column chosen by class name
    X = predictorExtractionFcn(tbl);
    [rawLabels, rawScores] = predict(model, X);
    % Find PVC column
    classNames = model.ClassNames; % {'Other','PVC'}
    if iscell(classNames)
        pvcIdx = find(strcmp(classNames, 'PVC'), 1);
    else
        pvcIdx = find(classNames == 'PVC', 1);
    end
    if isempty(pvcIdx)
        % Fallback: default to second column
        pvcIdx = min(size(rawScores,2), 2);
    end
    pvcScore = rawScores(:, pvcIdx);
    labels = cell(size(rawLabels));
    labels(:) = {'Other'};
    labels(pvcScore >= pvcThreshold) = {'PVC'};
    scores = rawScores;
end
