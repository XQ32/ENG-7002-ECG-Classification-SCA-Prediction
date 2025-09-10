function [ecg_filtered, filter_info] = ecgFilter(ecg, fs, method_index, power_line_freq)
% ========================================================================
% File: ecgFilter.m
% Overview: Unified entry for ECG preprocessing and multi-strategy baseline/noise suppression.
% Responsibilities:
%   1) Optional power-line notch (0=skip / 50 / 60 Hz) plus 2nd/3rd harmonic suppression.
%   2) Three baseline removal strategies: FIR high-pass / Butterworth high-pass /
%      enhanced pipeline (dual-median + high-pass).
%   3) Failure fallback chain: enhanced → IIR → FIR to ensure output under adverse conditions.
%   4) Low-pass to attenuate high-frequency and EMG noise (method 3 uses narrower bandwidth).
%   5) For method 3, additional morphology-preserving smoothing (SG / moving mean) and
%      QRS enhancement (differentiate-square-integrate).
%   6) Output complete filter parameters and intermediate flags for QA and tuning.
% Inputs:
%   ecg (double vector) : Raw single-lead ECG, row/column accepted, standardized to column.
%   fs (double)         : Sampling frequency in Hz.
%   method_index (int)  : 1=FIR high-pass, 2=Butterworth high-pass, 3=enhanced pipeline (default 1).
%   power_line_freq (int): Power-line frequency 0/50/60, 0 means skip (default 60).
% Outputs:
%   ecg_filtered (double column) : Final filtered signal (column vector).
%   filter_info (struct) : Method label, parameter struct, fallback flags, notch/lowpass/
%                          smoothing/QRS enhancement info.
% Key implementation notes:
%   - Power-line notch: narrowband IIR notch (Q≈35), add narrow 2nd/3rd harmonics if needed.
%   - FIR high-pass: kaiserord for adaptive order with cap ≤ N/10 to avoid overly high order.
%   - IIR high-pass: 4th-order Butterworth + filtfilt for zero-phase.
%   - Enhanced pipeline: dual-scale median baseline + 2nd-order high-pass to preserve morphology.
%   - QRS enhancement: classic Pan-Tompkins chain (diff→square→moving integrate).
%   - Consistent zero-phase processing to avoid ST/QRS distortion.
% Edge/robustness strategies:
%   - Invalid power-line value falls back to 60 Hz with a warning.
%   - Method 3/2 exceptions fall back stepwise (3→2→1), recording error / error2.
%   - Skip low-pass or harmonics when sampling rate is insufficient; use moving average when sgolay is unavailable.
%   - Windows and orders adapt to fs to reduce hard-coding.
% Changelog:
%   2025-08-30: Unified comment style rewrite; keep original algorithms/logic unchanged.
% ========================================================================

% Ensure ECG is a column vector for consistent processing
if size(ecg, 2) > size(ecg, 1)
    ecg = ecg';
end

% Default parameter settings
if nargin < 3
    method_index = 1;
end
if nargin < 4
    power_line_freq = 60; % default 60 Hz power-line interference
end

% Validate power-line parameter (0 allowed to skip notch)
if ~(power_line_freq == 0 || power_line_freq == 50 || power_line_freq == 60)
    warning('Invalid power-line frequency %d Hz; only 0(skip)/50/60 supported. Falling back to 60 Hz.', power_line_freq);
    power_line_freq = 60;
end

% High-pass cutoff set to 0.5 Hz
cutoff_freq = 0.5;

% Create output info struct
filter_info = struct();
filter_info.method_index = method_index;
filter_info.cutoff_freq = cutoff_freq;
filter_info.sampling_freq = fs;
filter_info.power_line_freq = power_line_freq;
filter_info.original = ecg;

%% Step 1: Apply notch filter to remove power-line interference (50/60 Hz)
if power_line_freq == 0
    % Skip notch (data already notch-filtered)
    ecg_notch_filtered = ecg;
    filter_info.notch_skipped = true;
    filter_info.has_2nd_harmonic_filter = false;
else
    % Design a narrowband power-line notch filter
    wo = power_line_freq/(fs/2);  % normalized notch frequency
    bw = wo/35;      % narrow bandwidth (Q factor = 35)
    [b, a] = iirnotch(wo, bw);  % IIR notch design

    % Apply notch with zero-phase filtering
    ecg_notch_filtered = filtfilt(b, a, ecg);
    filter_info.notch_filtered = ecg_notch_filtered;

    % Add power-line harmonic filters (2x / 3x)
    if fs > power_line_freq * 3 % ensure sampling rate high enough for harmonics
        % Add 2nd harmonic notch
        wo_2nd = (power_line_freq*2)/(fs/2);
        if wo_2nd < 0.95 % keep normalized frequency < 1
            bw_2nd = wo_2nd/40; % narrower harmonic bandwidth
            [b_2nd, a_2nd] = iirnotch(wo_2nd, bw_2nd);
            ecg_notch_filtered = filtfilt(b_2nd, a_2nd, ecg_notch_filtered);
            filter_info.has_2nd_harmonic_filter = true;
        end

        % Add 3rd harmonic notch
        wo_3rd = (power_line_freq*3)/(fs/2);
        if wo_3rd < 0.95
            bw_3rd = wo_3rd/45; % even narrower for 3rd harmonic
            [b_3rd, a_3rd] = iirnotch(wo_3rd, bw_3rd);
            ecg_notch_filtered = filtfilt(b_3rd, a_3rd, ecg_notch_filtered);
            filter_info.has_3rd_harmonic_filter = true;
        end
    end
end

%% Step 2: Choose baseline removal method by method_index
switch method_index
    case 1
        % Method 1: FIR high-pass filter
        [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
        filter_info.name = 'FIR high-pass filter';
        filter_info.filter_params = filter_params;
        
    case 2
        % Method 2: Butterworth IIR high-pass filter
        try
            [ecg_baseline_filtered, filter_params] = apply_IIR_filter(ecg_notch_filtered, fs, cutoff_freq);
            filter_info.name = 'Butterworth IIR high-pass filter';
            filter_info.filter_params = filter_params;
        catch e
            warning(e.identifier, 'Method 2 (IIR) failed: %s. Falling back to Method 1.', e.message);
            [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
            filter_info.name = 'FIR high-pass filter (fallback)';
            filter_info.filter_params = filter_params;
            filter_info.error = e.message;
        end
    
    case 3
        % Method 3: Enhanced denoise + smoothing pipeline (baseline removal + high-pass + strong low-pass)
        try
            [ecg_baseline_filtered, filter_params] = apply_enhanced_pipeline(ecg_notch_filtered, fs);
            filter_info.name = 'Enhanced denoise + smoothing pipeline';
            filter_info.filter_params = filter_params;
        catch e
            warning(e.identifier, 'Method 3 failed: %s. Falling back to Method 2.', e.message);
            try
                [ecg_baseline_filtered, filter_params] = apply_IIR_filter(ecg_notch_filtered, fs, cutoff_freq);
                filter_info.name = 'Butterworth IIR high-pass filter (fallback)';
                filter_info.filter_params = filter_params;
                filter_info.error = e.message;
            catch e2
                warning(e2.identifier, 'Method 2 also failed: %s. Finally falling back to Method 1.', e2.message);
                [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
                filter_info.name = 'FIR high-pass filter (final fallback)';
                filter_info.filter_params = filter_params;
                filter_info.error2 = e2.message;
            end
        end
        
    otherwise
        % Default to method 1
        warning('Invalid method index: %d. Using Method 1.', method_index);
        [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
        filter_info.name = 'FIR high-pass filter (default)';
        filter_info.filter_params = filter_params;
end

%% Step 3: Apply low-pass filter to remove high-frequency clutter
% Low-pass cutoff selection
% Methods 1/2 use wider bandwidth (150 Hz), method 3 uses stronger band limit (<= ~65 Hz)
if method_index == 3
    lowpass_cutoff = min(65, max(5, floor(fs/2 - 1))); % Hz
else
    lowpass_cutoff = 150; % Hz
end

% Design low-pass Butterworth filter
if fs > 2 * lowpass_cutoff % ensure Nyquist criterion
    lp_order = 4; % moderate roll-off
    [b_lp, a_lp] = butter(lp_order, lowpass_cutoff/(fs/2), 'low');
    ecg_filtered = filtfilt(b_lp, a_lp, ecg_baseline_filtered);
    
    % Record low-pass filter info
    filter_info.lowpass_applied = true;
    filter_info.lowpass_cutoff = lowpass_cutoff;
    filter_info.lowpass_params.order = lp_order;
    filter_info.lowpass_params.b = b_lp;
    filter_info.lowpass_params.a = a_lp;
    
    fprintf('  Applied %d Hz low-pass filter to remove high-frequency noise\n', lowpass_cutoff);
else
    % If sampling rate too low, skip low-pass
    ecg_filtered = ecg_baseline_filtered;
    filter_info.lowpass_applied = false;
    warning('Sampling rate %d Hz too low; skipping %d Hz low-pass filter', fs, lowpass_cutoff);
end

% Method 3: further morphology-preserving smoothing and QRS enhancement output
if method_index == 3
    try
        % Savitzky-Golay smoothing (preserve P/QRS/T morphology as much as possible)
        sg_order = 2;
        sg_frame = max(5, 2*floor(0.05*fs)+1); % ~80ms, 奇数
        if exist('sgolayfilt','file')
            ecg_smoothed = sgolayfilt(ecg_filtered, sg_order, sg_frame);
        else
            % Compatible fallback when sgolayfilt is unavailable: moving average
            ma_frame = max(3, round(0.05*fs));
            ecg_smoothed = movmean(ecg_filtered, ma_frame);
        end
        ecg_filtered = ecg_smoothed;
        filter_info.sgolay.order = sg_order;
        filter_info.sgolay.frame = sg_frame;
        filter_info.sgolay.applied = true;
    catch
        % Even if smoothing fails, do not affect main flow
        filter_info.sgolay.applied = false;
    end

    % Generate QRS-enhanced signal (Pan-Tompkins style: diff-square-moving integrate)
    try
        der_kernel = (1/8) * [1 2 0 -2 -1];
        qrs_diff = filtfilt(der_kernel, 1, ecg_filtered);
        qrs_squared = qrs_diff .^ 2;
        integ_win = max(5, round(0.12*fs));
        qrs_enhanced = movmean(qrs_squared, integ_win);
        filter_info.qrs_enhanced = qrs_enhanced;
        filter_info.qrs_enhance_params.diff_kernel = der_kernel;
        filter_info.qrs_enhance_params.window = integ_win;
    catch
        % If enhancement fails, provide empty field
        filter_info.qrs_enhanced = [];
    end
end

end

% Method 1: FIR high-pass filter function
function [filtered_signal, params] = apply_FIR_filter(signal, fs, cutoff_freq)
    % Get signal length
    signal_length = length(signal);
    
    % Compute acceptable max filter order (1/10 of signal length)
    max_filter_order = floor(signal_length/10);
    
    fcuts = [(cutoff_freq-0.07) (cutoff_freq)];
    mags = [0 1];
    % Relax tolerances to reduce required filter order
    devs = [0.01 0.005];
    
    % Design with kaiserord
    [n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fs);
    
    % Limit filter order to avoid overly high orders
    if n > max_filter_order
        % If calculated order too high, cap to lower order
        original_n = n;
        n = max_filter_order;
        fprintf('Filter order reduced from %d to %d\n', original_n, max_filter_order);
    end
    
    % Design FIR filter
    b = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');
    a = 1;
    
    % Apply zero-phase filtering
    filtered_signal = filtfilt(b, a, signal);
    
    % Record parameters
    params.filter_order = n;
    params.window = 'kaiser';
    params.beta = beta;
    params.cutoff = Wn;
    params.filter_type = ftype;
    params.b = b;
    params.a = a;
end

% Method 2: Butterworth IIR high-pass filter function
function [filtered_signal, params] = apply_IIR_filter(signal, fs, cutoff_freq)
    % Design Butterworth IIR high-pass filter
    order = 4;  % filter order
    [b, a] = butter(order, cutoff_freq/(fs/2), 'high');
    
    % Apply zero-phase filtering
    filtered_signal = filtfilt(b, a, signal);
    
    % Record parameters
    params.filter_order = order;
    params.filter_type = 'butterworth';
    params.cutoff = cutoff_freq;
    params.b = b;
    params.a = a;
end

% Method 3: Enhanced denoise + smoothing front pipeline (baseline removal + high-pass)
function [enhanced_signal, params] = apply_enhanced_pipeline(signal, fs)
    params = struct();

    % Dual-window median filtering to estimate and remove baseline (~0.2 s and ~0.6 s)
    win1 = max(3, 2*floor(0.20*fs)+1);
    win2 = max(3, 2*floor(0.60*fs)+1);
    if exist('medfilt1','file')
        baseline_rough = medfilt1(signal, win1, 'truncate');
        baseline = medfilt1(baseline_rough, win2, 'truncate');
    else
        % Compatible environment: fall back to moving median
        baseline_rough = movmedian(signal, win1);
        baseline = movmedian(baseline_rough, win2);
    end
    signal_detrended = signal - baseline;

    % Mild high-pass to further suppress drift (no forced band-pass; aligns with previous step)
    hp_cut = 0.5; % Hz
    [b_hp, a_hp] = butter(2, max(hp_cut, 0.01)/(fs/2), 'high');
    enhanced_signal = filtfilt(b_hp, a_hp, signal_detrended);

    % Record parameters
    params.baseline_median.win1 = win1;
    params.baseline_median.win2 = win2;
    params.highpass.order = 2;
    params.highpass.cutoff = hp_cut;
    params.highpass.b = b_hp;
    params.highpass.a = a_hp;
end
