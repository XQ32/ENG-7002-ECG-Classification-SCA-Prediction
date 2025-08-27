function [isGood, corrValues, usedTemplateClass, thresholdsOut] = assessBeatsQuality(segmentsCell, rIndices, beatTypes, fs, options)
% File: assessBeatsQuality.m
% Type: Function
% Usage:
%   [isGood, corrValues, usedTemplateClass, thresholdsOut] = assessBeatsQuality(segmentsCell, rIndices, beatTypes, fs, options)
%
% Description:
%   Use template correlation (SQI) to filter low-quality/abnormal beats. Build templates for
%   Other/PVC/Global, compute per-beat maximum correlation and compare with class-specific
%   thresholds, and output the pass mask and the template class used.
%
% Inputs:
%   - segmentsCell (cell{N}): per-beat segments (column vectors)
%   - rIndices (double[Nx1]): R index within each segment
%   - beatTypes (cellstr[Nx1]): beat type (e.g., 'Other'/'PVC')
%   - fs (double): sampling frequency (Hz)
%   - options (struct, optional): .windowSec/.bandpassHz/.minBeatsForTemplate/
%       .corrThreshold.(Other|PVC|Global)/.maxShiftSec
%
% Outputs:
%   - isGood (logical[Nx1]), corrValues (double[Nx1])
%   - usedTemplateClass (cellstr[Nx1])
%   - thresholdsOut (struct)
%
% Dependencies:
%   Local helpers: extractWindowAroundR, bandpassSafe, normalizeVec,
%                  maxCorrWithShift, mergeStruct
%
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26


    if nargin < 6 || isempty(options)
        options = struct();
    end

    % Default parameters
    defaults.windowSec = 0.24;
    defaults.bandpassHz = [5 40];
    defaults.minBeatsForTemplate = 12;
    defaults.corrThreshold = struct('Other', 0.80, 'PVC', 0.70, 'Global', 0.75);
    defaults.maxShiftSec = 0.010; % 10ms
    cfg = mergeStruct(defaults, options);

    N = numel(segmentsCell);
    isGood = true(N,1);
    corrValues = zeros(N,1);
    usedTemplateClass = repmat({''}, N, 1);
    thresholdsOut = cfg.corrThreshold;

    if N == 0
        return;
    end

    % Preprocessing: extract a fixed-length window centered at R, then band-pass + normalize
    L = max(4, 2*round(cfg.windowSec*fs/2)); % even length
    halfWin = L/2;
    maxShift = max(0, round(cfg.maxShiftSec * fs));

    X = zeros(L, N);
    validMask = true(N,1);
    for i = 1:N
        seg = segmentsCell{i};
        if isempty(seg) || ~isfinite(rIndices(i))
            validMask(i) = false; continue;
        end
        rIdx = rIndices(i);
        [w, ok] = extractWindowAroundR(seg, rIdx, halfWin);
        if ~ok
            validMask(i) = false; continue;
        end
        % Band-pass + normalization
        w = bandpassSafe(w, fs, cfg.bandpassHz);
        w = w - median(w);
        nrm = norm(w) + eps;
        w = w ./ nrm;
        X(:, i) = w(:);
    end

    if ~any(validMask)
        % All invalid → pass all (no gating)
        isGood = true(N,1);
        corrValues = zeros(N,1);
        usedTemplateClass(:) = {'Global'};
        return;
    end

    % Template classes: Other and PVC
    types = beatTypes(:);
    isPVC = ismember(lower(string(types)), lower("PVC"));
    isOther = ~isPVC; % Treat the rest as Other

    % Build templates (use median for robustness)
    templates = struct();
    templates.Other.available = false;
    templates.PVC.available = false;
    templates.Global.available = false;

    if sum(validMask & isOther) >= cfg.minBeatsForTemplate
        templates.Other.vec = median(X(:, validMask & isOther), 2);
        templates.Other.vec = normalizeVec(templates.Other.vec);
        templates.Other.available = true;
    end
    if sum(validMask & isPVC) >= cfg.minBeatsForTemplate
        templates.PVC.vec = median(X(:, validMask & isPVC), 2);
        templates.PVC.vec = normalizeVec(templates.PVC.vec);
        templates.PVC.available = true;
    end
    if sum(validMask) >= max(6, ceil(cfg.minBeatsForTemplate/2))
        templates.Global.vec = median(X(:, validMask), 2);
        templates.Global.vec = normalizeVec(templates.Global.vec);
        templates.Global.available = true;
    end

    % If no template available, pass all (avoid false rejections)
    if ~templates.Other.available && ~templates.PVC.available && ~templates.Global.available
        isGood = true(N,1);
        corrValues = zeros(N,1);
        usedTemplateClass(:) = {'Global'};
        return;
    end

    % Compute correlations and decide
    for i = 1:N
        if ~validMask(i)
            isGood(i) = true; % Not enough to evaluate, pass
            corrValues(i) = 0;
            usedTemplateClass{i} = 'Global';
            continue;
        end

        % Select template class
        if isPVC(i) && templates.PVC.available
            tmpl = templates.PVC.vec; tmplClass = 'PVC'; thr = cfg.corrThreshold.PVC;
        elseif ~isPVC(i) && templates.Other.available
            tmpl = templates.Other.vec; tmplClass = 'Other'; thr = cfg.corrThreshold.Other;
        elseif templates.Global.available
            tmpl = templates.Global.vec; tmplClass = 'Global'; thr = cfg.corrThreshold.Global;
        else
            % Should not reach here; allow by default
            isGood(i) = true; corrValues(i) = 0; usedTemplateClass{i} = 'Global';
            continue;
        end

        x = X(:, i);
        c = maxCorrWithShift(x, tmpl, maxShift);
        corrValues(i) = c;
        isGood(i) = (c >= thr);
        usedTemplateClass{i} = tmplClass;
    end
end

% ==== Helper functions ====
function [w, ok] = extractWindowAroundR(seg, rIdx, halfWin)
    L = 2*halfWin;
    n = numel(seg);
    s = round(rIdx - halfWin + 1);
    e = s + L - 1;
    ok = true;
    if s < 1 || e > n
        % Zero-pad at boundaries
        w = zeros(L,1);
        srcS = max(1, s);
        srcE = min(n, e);
        dstS = max(1, 1 + (srcS - s));
        dstE = dstS + (srcE - srcS);
        w(dstS:dstE) = seg(srcS:srcE);
    else
        w = seg(s:e);
    end
end

function v = bandpassSafe(x, fs, band)
    x = x(:);
    if numel(x) < 8
        v = x; return;
    end
    f1 = band(1); f2 = band(2);
    nyq = fs/2;
    if f2 >= nyq*0.95
        f2 = nyq*0.95;
    end
    if f1 <= 0
        f1 = 0.5;
    end
    if f2 <= f1
        v = x; return;
    end
    try
        [b, a] = butter(4, [f1, f2] / nyq, 'bandpass');
        v = filtfilt(b, a, x);
    catch
        v = x; % On failure, return input as-is
    end
end

function v = normalizeVec(v)
    v = v(:) - median(v);
    nrm = norm(v) + eps;
    v = v ./ nrm;
end

function c = maxCorrWithShift(x, y, maxShift)
    % x,y are already normalized
    if maxShift <= 0
        c = max(-1, min(1, x' * y));
        return;
    end
    c = -1;
    L = numel(x);
    for k = -maxShift:maxShift
        if k >= 0
            xs = [x(1+k:end); zeros(k,1)];
        else
            kk = -k;
            xs = [zeros(kk,1); x(1:end-kk)];
        end
        ck = xs' * y;
        if ck > c
            c = ck;
        end
    end
    c = max(-1, min(1, c));
end

function out = mergeStruct(base, opt)
    out = base;
    if ~isstruct(opt)
        return;
    end
    f = fieldnames(opt);
    for i = 1:numel(f)
        out.(f{i}) = opt.(f{i});
    end
end


