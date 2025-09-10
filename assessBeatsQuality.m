function [isGood, corrValues, usedTemplateClass, thresholdsOut] = assessBeatsQuality(segmentsCell, rIndices, beatTypes, fs, options)
% =========================================================================
% File name: assessBeatsQuality.m
% Overview: Build multi-class (QRS morphology) templates from beat segments
%           and compute a signal quality index (SQI) as the "max correlation
%           with small temporal shift". Filter beats by morphological quality,
%           produce a quality mask and correlation scores, and return the
%           template class and thresholds used.
%
% Responsibilities:
%   1. Preprocess each beat window (bandpass + median detrend + energy
%      normalization) and align length
%   2. Build "median templates" (robust centers) by class (Other / PVC) and
%      globally (Global)
%   3. For each beat, maximize correlation within ±maxShift samples via
%      circular-like shift
%   4. Apply class-wise or global correlation thresholds to decide quality
%      (pass/reject)
%   5. Use lenient pass-through when samples are insufficient / no valid
%      template exists, to avoid over-rejection
%
% Inputs:
%   segmentsCell : {N×1}  Each element is a column-vector beat segment
%                            (containing the R location)
%   rIndices     : N×1    R index within each segment (1-based)
%   beatTypes    : {N×1}  Beat type strings (e.g., 'PVC','Other')
%   fs           : double Sampling rate (Hz)
%   options      : struct Optional overrides:
%                  .windowSec             Total analysis window len (s,
%                                         centered at R)      default 0.24
%                  .bandpassHz            Bandpass [f1 f2] for preprocessing
%                                         default [5 40]
%                  .minBeatsForTemplate   Min beats to build a class template
%                                         default 12
%                  .corrThreshold.Other   Correlation threshold for Other
%                                         default 0.80
%                  .corrThreshold.PVC     Correlation threshold for PVC
%                                         default 0.70
%                  .corrThreshold.Global  Correlation threshold for Global
%                                         default 0.75
%                  .maxShiftSec           Allowed micro shift for maximizing
%                                         correlation (s)      default 0.010
%
% Outputs:
%   isGood            : N×1 logical  Quality pass flag (true=keep)
%   corrValues        : N×1 double   Max correlation with used template
%                                     under allowed shift
%   usedTemplateClass : {N×1}        The template class actually used
%                                     ('Other'/'PVC'/'Global')
%   thresholdsOut     : struct       The set of thresholds used (per class
%                                     and global)
%
% Dependencies (local helpers in this file):
%   extractWindowAroundR, bandpassSafe, normalizeVec,
%   maxCorrWithShift, mergeStruct
%
% Key implementation notes:
%   - Use median to build templates to suppress outliers
%   - Align length + unit-energy normalization to ensure fair comparison
%   - Allow ±maxShift shift to approximate the "cross-correlation peak"
%     via max dot product
%   - Multi-level fallback: insufficient class template → use Global;
%     none available → pass-through
%   - Safety: when data is empty or no valid template exists, do not force
%     rejection (avoid over-pruning)
%
% Changelog:
%   2025-08-30 Rewrite comments to unified style (no algorithm change)
% =========================================================================
%
% Input recap:
%   segmentsCell  - Cell array of beat segments (each a column vector)
%   rIndices      - R indices within segments (aligned with segmentsCell)
%   beatTypes     - Cell array of beat types ('Other' or 'PVC', etc.)
%   fs            - Sampling frequency (Hz)
%   options       - Optional config struct:
%                   .windowSec (default 0.24)   total window length centered at R
%                   .bandpassHz (default [5 40]) preprocessing bandpass
%                   .minBeatsForTemplate (default 12) min samples per class
%                   .corrThreshold.Other (default 0.80)
%                   .corrThreshold.PVC   (default 0.70)
%                   .corrThreshold.Global(default 0.75)
%                   .maxShiftSec (default 0.010) allowed micro alignment shift
%
% Output recap:
%   isGood            - logical vector, true=passes SQI
%   corrValues        - max correlation vs. the template used (with shift)
%   usedTemplateClass - template class used for each beat: 'Other'/'PVC'/'Global'
%   thresholdsOut     - thresholds struct used

    if nargin < 6 || isempty(options)
        options = struct();
    end

    % Defaults (tuned: improve true R recall, suppress T-as-R)
    defaults.windowSec = 0.30;                % was 0.24 → 0.30 s (±150 ms) to better cover wide QRS
    defaults.bandpassHz = [8 40];             % was [5 40] → [8 40] to suppress T/slow drift
    defaults.minBeatsForTemplate = 8;         % was 12 → 8; easier to form templates
    defaults.corrThreshold = struct( ...      % Correlation thresholds: slightly relax Global and Other
        'Other', 0.78, ...                    % was 0.80 → 0.78 (favor recall)
        'PVC',   0.70, ...                    % unchanged
        'Global',0.74  ...                    % was 0.75 → 0.74
    );
    defaults.maxShiftSec = 0.025; % was 10 ms → 25 ms; tolerate alignment errors, raise true-R correlation
    cfg = mergeStruct(defaults, options);

    N = numel(segmentsCell);
    isGood = true(N,1);
    corrValues = zeros(N,1);
    usedTemplateClass = repmat({''}, N, 1);
    thresholdsOut = cfg.corrThreshold;

    if N == 0
        return;
    end

    % Preprocess: extract fixed-length window centered at R; bandpass + normalize
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
    % Bandpass + normalization
        w = bandpassSafe(w, fs, cfg.bandpassHz);
        w = w - median(w);
        nrm = norm(w) + eps;
        w = w ./ nrm;
        X(:, i) = w(:);
    end

    if ~any(validMask)
        % All invalid → mark all as pass (no gating)
        isGood = true(N,1);
        corrValues = zeros(N,1);
        usedTemplateClass(:) = {'Global'};
        return;
    end

    % Class templates: Other and PVC
    types = beatTypes(:);
    isPVC = ismember(lower(string(types)), lower("PVC"));
    isOther = ~isPVC; % all others treated as Other

    % Build templates (median for robustness)
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

    % If no class/global templates are available, pass through (avoid over-rejection)
    if ~templates.Other.available && ~templates.PVC.available && ~templates.Global.available
        isGood = true(N,1);
        corrValues = zeros(N,1);
        usedTemplateClass(:) = {'Global'};
        return;
    end

    % Compute correlation and decide
    for i = 1:N
        if ~validMask(i)
            isGood(i) = true; % insufficient to assess; pass through
            corrValues(i) = 0;
            usedTemplateClass{i} = 'Global';
            continue;
        end

        % Choose template class
        if isPVC(i) && templates.PVC.available
            tmpl = templates.PVC.vec; tmplClass = 'PVC'; thr = cfg.corrThreshold.PVC;
        elseif ~isPVC(i) && templates.Other.available
            tmpl = templates.Other.vec; tmplClass = 'Other'; thr = cfg.corrThreshold.Other;
        elseif templates.Global.available
            tmpl = templates.Global.vec; tmplClass = 'Global'; thr = cfg.corrThreshold.Global;
        else
            % Should not reach here; fallback pass through
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
    % x,y are already unit-normalized
    if maxShift <= 0
        c = max(-1, min(1, x' * y));
        return;
    end
    c = -1;
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
