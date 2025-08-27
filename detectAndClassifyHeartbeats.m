function [heartbeatSegments, beatInfo, stats] = detectAndClassifyHeartbeats(ecg, ATRTIMED, ANNOTD, fs)
% File: detectAndClassifyHeartbeats.m
% Type: Function
% Usage:
%   [heartbeatSegments, beatInfo, stats] = detectAndClassifyHeartbeats(ecg, ATRTIMED, ANNOTD, fs)
%
% Description:
%   Robust R detection with mixed polarity/derivative energy/amplitude remapping, polarity-adaptive
%   Q/S localization, T-wave suppression, nearest-neighbor matching to R annotations, and SQI-based
%   quality filtering. Outputs segments and structured beat information.
%
% Inputs:
%   - ecg (double[Nx1]), ATRTIMED (double[Kx1]), ANNOTD (cellstr[Kx1]), fs (double)
%
% Outputs:
%   - heartbeatSegments (cell{M}), beatInfo (struct[M]), stats (struct)
%
% Dependencies:
%   assessBeatsQuality and local helper functions in this file
%
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26


    % Stats initialization
    stats = struct('num_detected_R', 0, ...
                   'num_r_matched_and_kept', 0, ...
                   'num_r_unmatched_removed', 0, ...
                   'num_matched_but_qs_removed', 0, ...
                   'mteo_failed', false, ...
                   'num_t_as_r_removed', 0);

    % Ensure column vector
    ecg = ecg(:);

    % Convert annotation times (s) to sample indices
    ann_sample_indices = round(ATRTIMED(:) * fs);
    ann_labels = ANNOTD(:);

    % Precompute (shared by R/Q/S): moving median baseline and PVC amplitude threshold
    base_win = max(5, round(1.2 * fs));
    try
        base_vec = movmedian(ecg, base_win);
    catch
        try
            base_vec = medfilt1(ecg, 1 + 2*round(0.6*fs));
        catch
            base_vec = zeros(size(ecg));
        end
    end
    abs_ecg_corrected = abs(ecg - base_vec);
    amp_pvc_thr = median(abs_ecg_corrected) + 3.5*local_mad(abs_ecg_corrected);
    dx_full = [0; diff(ecg)];
    dx_abs_full = abs(dx_full);
    thr_slope_full = median(dx_abs_full) + 2.5*local_mad(dx_abs_full);
    % Energy envelope (QRS gating) to avoid P/T being mistaken for Q/S
    der2 = dx_full.^2;
    env_win = max(1, round(0.06 * fs));
    energy_env = movmean(der2, env_win);
    env_thr = median(energy_env) + 1.5*local_mad(energy_env);

    % ========== Robust R-peak detection with mixed polarity + derivative energy + amplitude remapping (fast) ==========
    try
        r_indices = local_detect_r_hybrid_fast(ecg, fs, base_vec, dx_abs_full, thr_slope_full);
    catch ME
        fprintf('R-peak detection failed (%s), returning empty results.\n', ME.message);
        heartbeatSegments = cell(0,1);
        beatInfo = struct('beatType', {}, 'segment', {}, 'rIndex', {}, 'qIndex', {}, 'sIndex', {}, 'qrsOnIndex', {}, 'qrsOffIndex', {}, 'segmentStartIndex', {});
        stats.num_detected_R = 0;
        return;
    end

    if isempty(r_indices)
        heartbeatSegments = cell(0,1);
        beatInfo = struct('beatType', {}, 'segment', {}, 'rIndex', {}, 'qIndex', {}, 'sIndex', {}, 'qrsOnIndex', {}, 'qrsOffIndex', {}, 'segmentStartIndex', {});
        fprintf('No R-peaks detected, returning empty results.\n');
        stats.num_detected_R = 0;
        return;
    end

    % Keep indices within signal range
    r_indices = r_indices(r_indices >= 1 & r_indices <= length(ecg));
    stats.num_detected_R = numel(r_indices);

    % Post-hoc suppression: remove near-T peaks misdetected as R
    % Return T-wave "seed" indices (global), aligned with R, for later labeling
    [r_indices, t_seed_by_r, num_t_removed] = local_remove_t_as_r(ecg, r_indices, fs, base_vec, dx_full, energy_env);
    if ~isempty(num_t_removed)
        stats.num_t_as_r_removed = num_t_removed;
    end

    % If no annotations: keep all R peaks as beats, type set to 'Other'
    no_annotations = isempty(ann_sample_indices);

    % Prepare outputs (filled dynamically)
    kept_segments = cell(numel(r_indices), 1);
    kept_beats = repmat(struct('beatType', 'Other', 'segment', [], 'rIndex', 0, ...
                               'qIndex', NaN, 'sIndex', NaN, 'tIndex', NaN, 'qrsOnIndex', NaN, 'qrsOffIndex', NaN, 'segmentStartIndex', 0), numel(r_indices), 1);
    kept_count = 0;

    half_window = round(0.2 * fs);     % default left window 0.2 s (enough for Q)
    q_search_pre = round(0.20 * fs);   % search Q within 200 ms before R (PVC may have wider QRS)
    s_search_post = round(0.22 * fs);  % search S within 220 ms after R (PVC may have wider S)
    r_ann_match_tol = max(1, round(0.20 * fs));  % tolerance for nearest R-annotation match ±0.1 s

    % Precomputed: global slope threshold (shared) - already computed

    % Define loop upper bound
    num_beats = numel(r_indices);
    for i = 1:num_beats
        r_idx = r_indices(i);

        % Segment extraction: adaptively extend right window if PVC suspected
        base_here = base_vec(r_idx);
        r_val_here = ecg(r_idx) - base_here;
        s_win_ext_here = 0;
        if abs(r_val_here) > amp_pvc_thr
            s_win_ext_here = round(0.10 * fs);
        end
        half_window_right = max(half_window, s_search_post + s_win_ext_here + round(0.02*fs));
        seg_start = max(1, r_idx - half_window);
        seg_end   = min(length(ecg), r_idx + half_window_right);
        segment = ecg(seg_start:seg_end);
        r_in_seg = r_idx - seg_start + 1;

        % Nearest-neighbor matching to R annotations (if provided); otherwise keep
        nearest_ann_idx = [];
        if ~no_annotations
            [min_diff, nearest_ann_idx] = min(abs(ann_sample_indices - r_idx));
            if isempty(nearest_ann_idx) || isnan(min_diff) || min_diff > r_ann_match_tol
                stats.num_r_unmatched_removed = stats.num_r_unmatched_removed + 1;
                continue;
            end
        end

        % Q/S and QRS on/off: robust, polarity-adaptive local search (PVC-friendly)
        [q_idx_global, s_idx_global, qrs_on_global, qrs_off_global] = local_detect_qs(ecg, r_idx, fs, q_search_pre, s_search_post, base_vec, dx_abs_full, thr_slope_full, amp_pvc_thr, dx_full, energy_env, env_thr);

        % Map to segment indices
        if ~isnan(q_idx_global)
            q_in_seg = q_idx_global - seg_start + 1;
            if q_in_seg < 1 || q_in_seg > numel(segment)
                q_in_seg = NaN;
            end
        else
            q_in_seg = NaN;
        end

        if ~isnan(s_idx_global)
            s_in_seg = s_idx_global - seg_start + 1;
            if s_in_seg < 1 || s_in_seg > numel(segment)
                s_in_seg = NaN;
            end
        else
            s_in_seg = NaN;
        end

        % Map QRS on/off to segment
        if ~isnan(qrs_on_global)
            qrs_on_in_seg = qrs_on_global - seg_start + 1;
            if qrs_on_in_seg < 1 || qrs_on_in_seg > numel(segment)
                qrs_on_in_seg = NaN;
            end
        else
            qrs_on_in_seg = NaN;
        end
        if ~isnan(qrs_off_global)
            qrs_off_in_seg = qrs_off_global - seg_start + 1;
            if qrs_off_in_seg < 1 || qrs_off_in_seg > numel(segment)
                qrs_off_in_seg = NaN;
            end
        else
            qrs_off_in_seg = NaN;
        end

        % T wave detection (tag only, not used for features)
        t_seed_here = NaN;
        if exist('t_seed_by_r','var') && numel(t_seed_by_r) >= i
            t_seed_here = t_seed_by_r(i);
        end
        t_idx_global = local_detect_twave(ecg, r_idx, fs, base_vec, energy_env, env_thr, dx_full, qrs_off_global, t_seed_here);

        % Map T to segment
        if ~isnan(t_idx_global)
            t_in_seg = t_idx_global - seg_start + 1;
            if t_in_seg < 1 || t_in_seg > numel(segment)
                t_in_seg = NaN;
            end
        else
            t_in_seg = NaN;
        end

        % Determine type
        if no_annotations
            beat_type = 'Other';
        else
            beat_type = char(ann_labels{nearest_ann_idx});
        end

        % Fill outputs (keep)
        kept_count = kept_count + 1;
        kept_segments{kept_count,1} = segment;
        kept_beats(kept_count).beatType = beat_type;
        kept_beats(kept_count).segment = segment;
        kept_beats(kept_count).rIndex = r_in_seg;
        kept_beats(kept_count).qIndex = q_in_seg;
        kept_beats(kept_count).sIndex = s_in_seg;
        kept_beats(kept_count).tIndex = t_in_seg;
        kept_beats(kept_count).qrsOnIndex = qrs_on_in_seg;
        kept_beats(kept_count).qrsOffIndex = qrs_off_in_seg;
        kept_beats(kept_count).segmentStartIndex = seg_start;
        stats.num_r_matched_and_kept = stats.num_r_matched_and_kept + 1;
    end

    % Trim to candidates
    if kept_count > 0
        heartbeatSegments = kept_segments(1:kept_count);
        beatInfo = kept_beats(1:kept_count);
    else
        heartbeatSegments = cell(0,1);
        beatInfo = struct('beatType', {}, 'segment', {}, 'rIndex', {}, 'qIndex', {}, 'sIndex', {}, 'qrsOnIndex', {}, 'qrsOffIndex', {}, 'segmentStartIndex', {});
    end

    % ========== Template-correlation-based SQI filtering ==========
    % Purpose: quality assessment of candidate beats; remove unreliable/noisy beats
    num_candidates = numel(heartbeatSegments);
    if num_candidates > 0
        try
            % Prepare inputs
            segmentsCell = heartbeatSegments;
            rIndices = arrayfun(@(b) double(b.rIndex), beatInfo);
            beatTypes = arrayfun(@(b) char(string(b.beatType)), beatInfo, 'UniformOutput', false);

            % Run SQI (default parameters)
            [isGood, corrValues, usedTemplateClass, thresholdsOut] = assessBeatsQuality(segmentsCell, rIndices, beatTypes, fs, struct()); %#ok<ASGLU>

            % Write SQI back
            for ii = 1:num_candidates
                beatInfo(ii).sqiIsGood = logical(isGood(ii));
                beatInfo(ii).sqiCorr = corrValues(ii);
                beatInfo(ii).sqiTemplateClass = usedTemplateClass{ii};
            end

            % Filter by SQI
            mask_keep = logical(isGood(:));
            num_removed_sqi = sum(~mask_keep);
            if any(~mask_keep)
                heartbeatSegments = heartbeatSegments(mask_keep);
                beatInfo = beatInfo(mask_keep);
            end

            % Stats
            if ~isfield(stats, 'num_sqi_removed')
                stats.num_sqi_removed = 0;
            end
            stats.num_sqi_removed = stats.num_sqi_removed + num_removed_sqi;
            if isfield(stats, 'num_r_matched_and_kept')
                stats.num_r_matched_and_kept = numel(beatInfo);
            end
        catch ME
            fprintf('SQI evaluation failed (%s), skipping quality gating and keeping all candidate beats.\n', ME.message);
        end
    end

    % Output stats
    if ~isempty(beatInfo)
        fprintf('Processing complete. Kept %d beats (R±0.2s).\n', numel(beatInfo));
    else
        fprintf('Processing complete. No qualifying beats were retained.\n');
    end

end

% ======================== Local helper functions ========================

function r_indices = local_detect_r_hybrid_fast(x, fs, base_vec, dx_abs_full, thr_slope_full)
% Fast R-peak detection: abs amplitude coarse screening + derivative energy + |x| remapping;
% T-wave suppression; fallback when counts look abnormal
    x = double(x(:));
    n = numel(x);
    if n < max(10, round(0.5*fs))
        r_indices = [];
        return;
    end

    refractory = max(1, round(0.25 * fs));   % 250 ms minimum RR interval
    search_half = max(1, round(0.10 * fs));  % remapping window ±100 ms

    % 0) Coarse screening by absolute amplitude
    if nargin < 3 || isempty(base_vec), base_vec = zeros(size(x)); end
    absx_corrected = abs(x - base_vec);
    thr_abs0 = median(absx_corrected) + 2.5*local_mad(absx_corrected);
    [~, loc_abs0] = findpeaks(absx_corrected, 'MinPeakDistance', refractory, 'MinPeakHeight', thr_abs0);

    % 1) Squared derivative with moving integration (fast)
    dx = diff(x); dx2 = dx.^2;
    win = max(1, round(0.10 * fs));
    integ = movmean(dx2, win);
    thr_integ = median(integ) + 3*local_mad(integ);
    [~, loc_i1] = findpeaks(integ, 'MinPeakDistance', refractory, 'MinPeakHeight', thr_integ);

    % 2) Remap to |x| maxima
    r_from_i1 = local_refine_to_abs_peak(x, loc_i1, search_half);

    % Initial set: coarse abs ∪ remapped energy peaks
    cand = sort(unique([loc_abs0(:); r_from_i1(:)]));
    cand = cand(cand >= 1 & cand <= n);

    % Sanity check (too few/many) → add fallback candidates
    durSec = n / fs;
    if isempty(cand) || numel(cand) < max(1, floor(0.5 * durSec)) || numel(cand) > ceil(3.5 * durSec)
        absx = abs(x);
        thr_abs = median(absx) + 3*local_mad(absx);
        [~, loc_abs] = findpeaks(absx, 'MinPeakDistance', refractory, 'MinPeakHeight', thr_abs);
        cand = sort(unique([cand(:); loc_abs(:)]));

        % Still too few → try Pan-Tompkins as final fallback
        if numel(cand) < max(1, floor(0.4 * durSec))
            try
                [~, r_pt] = heplab_pan_tompkin(x, fs, 0);
                cand = sort(unique([cand(:); r_pt(:)]));
            catch
            end
        end
    end

    if isempty(cand)
        r_indices = cand; return;
    end

    % Merge close neighbors (keep larger |x|)
    r_indices = local_merge_close_by_amplitude(x, cand, refractory);

    % Light screening: slope threshold and morphology ratio (avoid T-waves)
    if numel(r_indices) > 1
        if nargin < 4 || isempty(dx_abs_full)
            dx_abs = [abs(diff(x)); 0];
        else
            dx_abs = dx_abs_full;
        end
        if nargin < 5 || isempty(thr_slope_full)
            thr_slope = median(dx_abs) + 3*local_mad(dx_abs);
        else
            thr_slope = max(thr_slope_full, median(dx_abs) + 2.5*local_mad(dx_abs));
        end
        keep = false(size(r_indices));
        slope_win = max(1, round(0.05*fs));
        for ii = 1:numel(r_indices)
            a = max(1, r_indices(ii)-slope_win);
            b = min(n, r_indices(ii)+slope_win);
            % Slope must be steep enough
            cond_slope = max(dx_abs(a:b)) >= thr_slope;
            % Morphology ratio within [R-0.08s,R] or [R,R+0.08s]
            w = max(1, round(0.08*fs));
            left = max(1, r_indices(ii)-w);
            right = min(n, r_indices(ii)+w);
            loc_amp = absx_corrected(r_indices(ii));
            neigh_med = median(absx_corrected([left:r_indices(ii)-1, min(n,r_indices(ii)+1):right]));
            if isempty(neigh_med) || isnan(neigh_med) || neigh_med==0
                neigh_med = median(absx_corrected(max(1,r_indices(ii)-2*w):min(n,r_indices(ii)+2*w)));
            end
            cond_ratio = loc_amp >= 1.6 * neigh_med; % suppress smooth T-wave peaks
            keep(ii) = cond_slope && cond_ratio;
        end
        if any(~keep)
            r_indices = r_indices(keep);
        end
    end
end

function r_ref = local_refine_to_abs_peak(x, locs, halfw)
% Map candidates within locs±halfw to the maximum of |x|
    x = x(:); n = numel(x);
    if isempty(locs), r_ref = locs; return; end
    if nargin < 3, halfw = 0; end
    r_ref = zeros(size(locs));
    for k = 1:numel(locs)
        a = max(1, locs(k)-halfw);
        b = min(n, locs(k)+halfw);
        [~, rel] = max(abs(x(a:b)));
        r_ref(k) = a + rel - 1;
    end
    r_ref = r_ref(:);
end

function out = local_merge_close_by_amplitude(x, idx, minDist)
% Merge too-close peaks: within minDist/2 keep the larger |x|
    if isempty(idx), out = idx; return; end
    idx = sort(idx(:));
    keep = true(size(idx));
    for i = 2:numel(idx)
        if (idx(i) - idx(i-1)) < floor(minDist/2)
            % Two are too close; keep the one with larger amplitude
            if abs(x(idx(i))) >= abs(x(idx(i-1)))
                keep(i-1) = false;
            else
                keep(i) = false;
            end
        end
    end
    out = idx(keep);
end

function m = local_mad(v)
% Simple MAD (independent of statistics toolbox)
    v = v(:);
    med = median(v);
    m = median(abs(v - med)) + eps;
end

function [q_idx, s_idx, qrs_on, qrs_off] = local_detect_qs(x, r_idx, fs, q_pre_win, s_post_win, base_vec, dx_abs_full, thr_slope, amp_pvc_thr, dx_full, energy_env, env_thr)
% Polarity-adaptive Q/S detection; search for extrema opposite to R polarity around R;
% allow wider window for PVC; estimate QRS on/off based on derivative thresholds.
    n = numel(x);
    r_idx = max(1, min(n, r_idx));
    if nargin < 4 || isempty(q_pre_win), q_pre_win = round(0.20*fs); end
    if nargin < 5 || isempty(s_post_win), s_post_win = round(0.22*fs); end
    if nargin < 6 || isempty(base_vec), base_vec = zeros(size(x)); end
    if nargin < 7 || isempty(dx_abs_full), dx_abs_full = [abs(diff(x)); 0]; end
    if nargin < 8 || isempty(thr_slope)
        thr_slope = median(dx_abs_full) + 2.5*local_mad(dx_abs_full);
    end
    if nargin < 9 || isempty(amp_pvc_thr)
        amp_pvc_thr = median(abs(x-base_vec)) + 3.5*local_mad(abs(x-base_vec));
    end
    if nargin < 10 || isempty(dx_full), dx_full = [0; diff(x)]; end
    if nargin < 11 || isempty(energy_env)
        der2 = dx_full.^2; energy_env = movmean(der2, max(1, round(0.06*fs)));
    end
    if nargin < 12 || isempty(env_thr)
        env_thr = median(energy_env) + 2.0*local_mad(energy_env);
    end

    % Use precomputed moving median baseline
    base = base_vec(r_idx);
    r_val = x(r_idx) - base;
    r_sign = sign(r_val); if r_sign == 0, r_sign = 1; end

    % --- Q search (before R; energy gate + fallback) ---
    pre_a = max(1, r_idx - q_pre_win);
    pre_b = max(1, r_idx - 1);
    q_idx = NaN;
    if pre_b >= pre_a
        gate = energy_env(pre_a:pre_b) >= env_thr * 0.8; % Q has lower energy; slightly relaxed threshold
        idxs = find(gate);
        if ~isempty(idxs)
            seg_all = x(pre_a:pre_b) - base;
            seg = seg_all(idxs);
            if r_sign > 0
                [~, rel] = min(seg);
            else
                [~, rel] = max(seg);
            end
            q_idx = pre_a + idxs(rel) - 1;
        else
            % Fallback: search the whole window
            seg = x(pre_a:pre_b) - base;
            if r_sign > 0, [~, rel] = min(seg); else, [~, rel] = max(seg); end
            q_idx = pre_a + rel - 1;
        end
        % If Q amplitude is too small, fallback to steepest slope point with opposite polarity in [R-0.16s, R-0.02s]
        if ~isnan(q_idx)
            if abs(x(q_idx)-base) < 0.08*abs(r_val)
                a2 = max(1, r_idx - round(0.16*fs));
                b2 = max(1, r_idx - round(0.02*fs));
                if b2 > a2
                    segdx = abs(dx_full(a2:b2));
                    [~, kmax] = max(segdx);
                    cand = a2 + kmax - 1;
                    if r_sign*(x(cand)-base) < 0
                        q_idx = cand;
                    end
                end
            end
        end
    end

    % --- S search (after R; energy gate + QRS off prior + PVC window extension + zero-cross/slope fallback + T protection) ---
    s_win_ext = 0;
    if abs(r_val) > amp_pvc_thr
        s_win_ext = round(0.10*fs); % PVC suspicion → extend S window
    end
    post_a = min(n, r_idx + 1);
    post_b = min(n, r_idx + s_post_win + s_win_ext);

    % Reference energy at R
    r_env = energy_env(r_idx);

    % Estimate QRS off using energy envelope, limit S search not to go beyond QRS
    qrs_off_env = NaN;
    max_search_env = min(n, r_idx + round(0.40*fs) + s_win_ext);
    dyn_thr_low = max(env_thr * 0.8, 0.25 * r_env);
    min_low_dur = max(1, round(0.04 * fs));
    cnt = 0;
    for tt = r_idx+1:max_search_env
        if energy_env(tt) < dyn_thr_low
            cnt = cnt + 1;
            if cnt >= min_low_dur
                qrs_off_env = tt;
                break;
            end
        else
            cnt = 0;
        end
    end
    if ~isnan(qrs_off_env)
        post_b = min(post_b, qrs_off_env + max(1, round(0.02*fs)));
    end
    s_idx = NaN;
    if post_b >= post_a
        % Dynamic energy gate using R energy as reference
        r_env = energy_env(r_idx);
        dyn_thr = max(env_thr, 0.35 * r_env);
        gate = energy_env(post_a:post_b) >= dyn_thr; % S likely in higher energy band
        idxs = find(gate);
        if ~isempty(idxs)
            seg_all = x(post_a:post_b) - base;
            seg = seg_all(idxs);
            if r_sign > 0, [~, rel] = min(seg); else, [~, rel] = max(seg); end
            s_idx = post_a + idxs(rel) - 1;
        else
            % Fallback: search whole window
            seg = x(post_a:post_b) - base;
            if r_sign > 0, [~, rel] = min(seg); else, [~, rel] = max(seg); end
            s_idx = post_a + rel - 1;
        end
        % Fallback: if same polarity as R or too small amplitude → find first opposite-polarity extremum after derivative zero-crossing
        if ~isnan(s_idx)
            if r_sign * (x(s_idx)-base) > 0 || abs(x(s_idx)-base) < 0.1*abs(r_val)
                % Find derivative zero-crossing
                search_a = r_idx;
                search_b = min(n, r_idx + round(0.30*fs) + s_win_ext);
                zc = find(dx_full(search_a:search_b-1).*dx_full(search_a+1:search_b) <= 0, 1, 'first');
                if ~isempty(zc)
                    zpos = search_a + zc - 1;
                    seg2 = x(zpos:search_b) - base;
                    if r_sign > 0
                        [~, rel2] = min(seg2);
                    else
                        [~, rel2] = max(seg2);
                    end
                    s_idx = zpos + rel2 - 1;
                else
                    % If no zero-crossing, fallback to opposite-polarity extremum near max slope
                    segdx2 = abs(dx_full(search_a:search_b));
                    [~, kmax2] = max(segdx2);
                    c2 = search_a + kmax2 - 1;
                    win = max(1, round(0.08*fs));
                    aa = max(1, c2 - win); bb = min(n, c2 + win);
                    seg3 = x(aa:bb) - base;
                    if r_sign > 0
                        [~, rr] = min(seg3);
                    else
                        [~, rr] = max(seg3);
                    end
                    s_idx = aa + rr - 1;
                end
            end
            % T-protection: if chosen S is too late relative to R, shorten to previous valley
            max_s_latency = round(0.27 * fs);
            if ~isnan(s_idx) && (s_idx - r_idx) > max_s_latency
                cut_b = r_idx + max_s_latency;
                seg4 = x(post_a:cut_b) - base;
                if r_sign > 0, [~, rr2] = min(seg4); else, [~, rr2] = max(seg4); end
                s_idx = post_a + rr2 - 1;
            end
        end
    end

    % --- QRS on/off by derivative threshold ---
    dx_abs = dx_abs_full;

    % Onset: from Q leftwards, last sample where slope stays below threshold
    qrs_on = NaN;
    if ~isnan(q_idx)
        left_a = max(1, q_idx - round(0.25*fs));
        low_mask = dx_abs(left_a:q_idx) < thr_slope;
        last_low = find(low_mask, 1, 'last');
        if ~isempty(last_low)
            qrs_on = left_a + last_low - 1;
        else
            qrs_on = left_a;
        end
    end

    % Offset: from S rightwards, first sample where slope stays below threshold
    qrs_off = NaN;
    if ~isnan(s_idx)
        right_b = min(n, s_idx + round(0.30*fs));
        low_mask = dx_abs(s_idx:right_b) < thr_slope;
        first_low = find(low_mask, 1, 'first');
        if ~isempty(first_low)
            qrs_off = s_idx + first_low - 1;
        else
            qrs_off = right_b;
        end
    end

    % Consistency checks
    if ~isnan(qrs_on) && ~isnan(q_idx)
        qrs_on = min(qrs_on, q_idx);
    end
    if ~isnan(qrs_off) && ~isnan(s_idx)
        qrs_off = max(qrs_off, s_idx);
    end
    if ~isnan(q_idx) && q_idx > r_idx, q_idx = NaN; end
    if ~isnan(s_idx) && s_idx < r_idx, s_idx = NaN; end
end

function [r_kept, t_seed_kept, num_removed] = local_remove_t_as_r(x, r_idx_list, fs, base_vec, dx_full, energy_env)
% Identify and remove R candidates that are actually T-waves.
% Heuristic: peaks 120–460 ms after R with clearly smaller energy/slope than the previous R are likely T-waves.
    tmp_sorted = sort(r_idx_list(:));
    r_idx_list = tmp_sorted.';
    n = numel(r_idx_list);
    if n <= 1
        r_kept = r_idx_list; t_seed_kept = NaN(size(r_kept)); num_removed = 0; return; end

    slope_win = max(1, round(0.06*fs));
    min_dt = round(0.12*fs); max_dt = round(0.46*fs);
    keep = true(1, n);
    t_seed = NaN(1, n);

    for k = 1:n-1
        r0 = r_idx_list(k); r1 = r_idx_list(k+1);
        dt = r1 - r0;
        if dt < min_dt || dt > max_dt
            continue; % outside T-wave window
        end
        base0 = base_vec(r0);
        r0_amp = abs(x(r0) - base0);
        r1_amp = abs(x(r1) - base_vec(r1));
        % Local slope
        a0 = max(1, r0 - slope_win); b0 = min(numel(x), r0 + slope_win);
        a1 = max(1, r1 - slope_win); b1 = min(numel(x), r1 + slope_win);
        s0 = max(abs(dx_full(a0:b0)));
        s1 = max(abs(dx_full(a1:b1)));
        % Energy envelope (vs R)
        e0 = energy_env(r0) + eps; e1 = energy_env(r1);
        ratio_amp = r1_amp / max(r0_amp, eps);
        ratio_slope = s1 / max(s0, eps);
        ratio_env = e1 / e0;
        % Same polarity looks more like T (not hard requirement)
        same_polar = sign(x(r0)-base0) * sign(x(r1)-base_vec(r1));
        if same_polar == 0, same_polar = 1; end
        % Decision: energy and slope both much smaller, and amplitude not too large
        if (ratio_env < 0.60 && ratio_slope < 0.65 && ratio_amp < 1.05 && same_polar > 0)
            keep(k+1) = false;      % treat the (k+1)-th as T-wave → remove
            t_seed(k) = r1;         % set as T-wave seed for the previous R
        end
    end

    r_kept = r_idx_list(keep);
    % T seeds aligned with the kept R list
    t_seed_kept = t_seed(keep);
    num_removed = sum(~keep);
end

function t_idx = local_detect_twave(x, r_idx, fs, base_vec, energy_env, env_thr, dx_full, qrs_off_global, seed)
% Detect T-wave peaks (tag only).
% Search window: max(R+80ms, QRSoff+20ms) to R+480ms;
% choose the same polarity as R; verify with slope/energy vs R.
    n = numel(x);
    r_idx = max(1, min(n, r_idx));
    if nargin < 4 || isempty(base_vec), base_vec = zeros(size(x)); end
    if nargin < 5 || isempty(energy_env)
        dx = [0; diff(x)]; der2 = dx.^2; energy_env = movmean(der2, max(1, round(0.06*fs)));
    end
    if nargin < 6 || isempty(env_thr)
        env_thr = median(energy_env) + 1.5*local_mad(energy_env);
    end
    if nargin < 7 || isempty(dx_full), dx_full = [0; diff(x)]; end

    base = base_vec(r_idx);
    r_val = x(r_idx) - base; r_sign = sign(r_val); if r_sign == 0, r_sign = 1; end
    start_idx = r_idx + round(0.08*fs);
    if nargin >= 8 && ~isnan(qrs_off_global)
        start_idx = max(start_idx, min(n, qrs_off_global + round(0.02*fs)));
    end
    stop_idx = min(n, r_idx + round(0.48*fs));
    if start_idx >= stop_idx
        t_idx = NaN; return; end

    % Prefer seed
    cand = NaN;
    if nargin >= 9 && ~isnan(seed)
        if seed >= start_idx- round(0.02*fs) && seed <= stop_idx + round(0.02*fs)
            cand = seed;
        end
    end
    if isnan(cand)
        seg = x(start_idx:stop_idx) - base;
        if r_sign > 0
            [~, rel] = max(seg);
        else
            [~, rel] = min(seg);
        end
        cand = start_idx + rel - 1;
    end

    % Verify: slope/energy should be clearly smaller than R
    win = max(1, round(0.06*fs));
    a0 = max(1, r_idx - win); b0 = min(n, r_idx + win);
    a1 = max(1, cand - win); b1 = min(n, cand + win);
    s0 = max(abs(dx_full(a0:b0))); s1 = max(abs(dx_full(a1:b1)));
    e0 = energy_env(r_idx) + eps; e1 = energy_env(cand);
    amp_ratio = abs(x(cand)-base) / max(abs(r_val), eps);
    if s1 < 0.70*s0 && e1 < 0.65*e0 && amp_ratio < 1.10
        t_idx = cand;
    else
        % Relax once but require low energy to avoid false positives
        gate = energy_env(start_idx:stop_idx) < max(env_thr, 0.5*e0);
        if isempty(find(gate, 1))
            t_idx = NaN; return; end
        seg2 = x(start_idx:stop_idx) - base;
        seg2(~gate) = -inf * sign(r_sign); % avoid high-energy points
        if r_sign > 0
            [~, rel2] = max(seg2);
        else
            [~, rel2] = min(seg2);
        end
        cand2 = start_idx + rel2 - 1;
        a1 = max(1, cand2 - win); b1 = min(n, cand2 + win);
        s1 = max(abs(dx_full(a1:b1))); e1 = energy_env(cand2);
        amp_ratio = abs(x(cand2)-base) / max(abs(r_val), eps);
        if s1 < 0.75*s0 && e1 < 0.70*e0 && amp_ratio < 1.10
            t_idx = cand2;
        else
            t_idx = NaN;
        end
    end
end