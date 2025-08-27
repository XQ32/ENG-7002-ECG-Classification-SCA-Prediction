function [ecg_filtered, filter_info] = ecgFilter(ecg, fs, method_index, power_line_freq)
% File: ecgFilter.m
% Type: Function
% Usage:
%   [ecg_filtered, filter_info] = ecgFilter(ecg, fs, method_index, power_line_freq)
%
% Description:
%   ECG preprocessing pipeline: optional mains notch (0/50/60 Hz) + baseline wander removal
%   (FIR/IIR/enhanced pipeline) + low-pass denoising. In method 3, optional morphology-preserving
%   smoothing and QRS enhancement are applied.
%
% Inputs:
%   - ecg (double[Vec]): raw ECG (row or column vector)
%   - fs (double): sampling rate (Hz)
%   - method_index (int): 1=FIR high-pass; 2=Butterworth IIR; 3=enhanced pipeline
%   - power_line_freq (int): 0(skip)/50/60 (default 60)
%
% Outputs:
%   - ecg_filtered (double[Col]): filtered ECG (column vector)
%   - filter_info (struct): filter configuration and parameters
%
% Dependencies:
%   Local helpers: apply_FIR_filter, apply_IIR_filter, apply_enhanced_pipeline
%
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

% Ensure ECG is a column vector for consistent processing
if size(ecg, 2) > size(ecg, 1)
    ecg = ecg';
end

% Default parameters
if nargin < 3
    method_index = 1;
end
if nargin < 4
    power_line_freq = 60; % default 60 Hz mains
end

% Validate mains parameter (allow 0 to skip notch)
if ~(power_line_freq == 0 || power_line_freq == 50 || power_line_freq == 60)
    warning('Invalid mains freq %d Hz. Only 0/50/60 supported. Using 60 Hz.', power_line_freq);
    power_line_freq = 60;
end

% Cutoff freq for baseline removal = 0.5 Hz
cutoff_freq = 0.5;

% Create output info struct
filter_info = struct();
filter_info.method_index = method_index;
filter_info.cutoff_freq = cutoff_freq;
filter_info.sampling_freq = fs;
filter_info.power_line_freq = power_line_freq;
filter_info.original = ecg;

%% Step 1: Notch filter for mains (50/60 Hz)
if power_line_freq == 0
    % Skip notch (data already notched at 60 Hz)
    ecg_notch_filtered = ecg;
    filter_info.notch_skipped = true;
    filter_info.has_2nd_harmonic_filter = false;
else
    % Design a narrowband mains notch filter
    wo = power_line_freq/(fs/2);  % normalized freq
    bw = wo/35;      % narrow bandwidth (Q=35)
    [b, a] = iirnotch(wo, bw);

    % Zero-phase filtering for notch
    ecg_notch_filtered = filtfilt(b, a, ecg);
    filter_info.notch_filtered = ecg_notch_filtered;

    % Add harmonic notches (2nd/3rd) if feasible
    if fs > power_line_freq * 3 % ensure Nyquist margin
        % 2nd harmonic
        wo_2nd = (power_line_freq*2)/(fs/2);
        if wo_2nd < 0.95 % ensure normalized freq < 1
            bw_2nd = wo_2nd/40; % narrower for harmonic
            [b_2nd, a_2nd] = iirnotch(wo_2nd, bw_2nd);
            ecg_notch_filtered = filtfilt(b_2nd, a_2nd, ecg_notch_filtered);
            filter_info.has_2nd_harmonic_filter = true;
        end

        % 3rd harmonic
        wo_3rd = (power_line_freq*3)/(fs/2);
        if wo_3rd < 0.95
            bw_3rd = wo_3rd/45; % even narrower
            [b_3rd, a_3rd] = iirnotch(wo_3rd, bw_3rd);
            ecg_notch_filtered = filtfilt(b_3rd, a_3rd, ecg_notch_filtered);
            filter_info.has_3rd_harmonic_filter = true;
        end
    end
end

%% Step 2: Baseline removal method by method_index
switch method_index
    case 1
        % Method 1: FIR high-pass
        [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
        filter_info.name = 'FIR High-pass';
        filter_info.filter_params = filter_params;
        
    case 2
        % Method 2: Butterworth IIR
        try
            [ecg_baseline_filtered, filter_params] = apply_IIR_filter(ecg_notch_filtered, fs, cutoff_freq);
            filter_info.name = 'Butterworth IIR';
            filter_info.filter_params = filter_params;
        catch e
            warning(e.identifier, 'Method 2 (IIR) failed: %s, fallback to method 1', e.message);
            [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
            filter_info.name = 'FIR High-pass (fallback)';
            filter_info.filter_params = filter_params;
            filter_info.error = e.message;
        end
    
    case 3
        % Method 3: Enhanced denoise + smoothing pipeline (baseline removal + HP + LP)
        try
            [ecg_baseline_filtered, filter_params] = apply_enhanced_pipeline(ecg_notch_filtered, fs);
            filter_info.name = 'Enhanced denoise + smoothing';
            filter_info.filter_params = filter_params;
        catch e
            warning(e.identifier, 'Method 3 failed: %s, fallback to method 2', e.message);
            try
                [ecg_baseline_filtered, filter_params] = apply_IIR_filter(ecg_notch_filtered, fs, cutoff_freq);
                filter_info.name = 'Butterworth IIR (fallback)';
                filter_info.filter_params = filter_params;
                filter_info.error = e.message;
            catch e2
                warning(e2.identifier, 'Method 2 also failed: %s, final fallback to method 1', e2.message);
                [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
                filter_info.name = 'FIR High-pass (final fallback)';
                filter_info.filter_params = filter_params;
                filter_info.error2 = e2.message;
            end
        end
        
    otherwise
        % Default to method 1
        warning('Invalid method index: %d, using method 1', method_index);
        [ecg_baseline_filtered, filter_params] = apply_FIR_filter(ecg_notch_filtered, fs, cutoff_freq);
        filter_info.name = 'FIR High-pass (default)';
        filter_info.filter_params = filter_params;
end

%% Step 3: Low-pass filter to remove high-frequency noise
% Low-pass cutoff selection:
% method 1/2: wider band (150 Hz); method 3: stronger limit (<=65 Hz)
if method_index == 3
    lowpass_cutoff = min(65, max(5, floor(fs/2 - 1))); % Hz
else
    lowpass_cutoff = 150; % Hz
end

% Design Butterworth low-pass
if fs > 2 * lowpass_cutoff % ensure Nyquist criterion
    lp_order = 4; % moderate roll-off
    [b_lp, a_lp] = butter(lp_order, lowpass_cutoff/(fs/2), 'low');
    ecg_filtered = filtfilt(b_lp, a_lp, ecg_baseline_filtered);
    
    % Record low-pass info
    filter_info.lowpass_applied = true;
    filter_info.lowpass_cutoff = lowpass_cutoff;
    filter_info.lowpass_params.order = lp_order;
    filter_info.lowpass_params.b = b_lp;
    filter_info.lowpass_params.a = a_lp;
    
    fprintf('  Applied %d Hz low-pass to remove high-frequency noise\n', lowpass_cutoff);
else
    % If sampling rate too low, skip low-pass
    ecg_filtered = ecg_baseline_filtered;
    filter_info.lowpass_applied = false;
    warning('Sampling rate %d Hz too low, skip %d Hz low-pass', fs, lowpass_cutoff);
end

% Method 3: further morphology-preserving smoothing and QRS enhancement
if method_index == 3
    try
        % Savitzky-Golay smoothing (preserve P/QRS/T morphology)
        sg_order = 2;
        sg_frame = max(5, 2*floor(0.05*fs)+1); % ~80ms, odd
        if exist('sgolayfilt','file')
            ecg_smoothed = sgolayfilt(ecg_filtered, sg_order, sg_frame);
        else
            % Without sgolayfilt, fallback to moving average
            ma_frame = max(3, round(0.05*fs));
            ecg_smoothed = movmean(ecg_filtered, ma_frame);
        end
        ecg_filtered = ecg_smoothed;
        filter_info.sgolay.order = sg_order;
        filter_info.sgolay.frame = sg_frame;
        filter_info.sgolay.applied = true;
    catch
        % Even if smoothing fails, main flow continues
        filter_info.sgolay.applied = false;
    end

    % Generate QRS-enhanced signal (Pan-Tompkins style: diff-square-integrate)
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
    % Signal length
    signal_length = length(signal);
    
    % Max acceptable filter order (1/10 of signal length)
    max_filter_order = floor(signal_length/10);
    
    fcuts = [(cutoff_freq-0.07) (cutoff_freq)];
    mags = [0 1];
    % Relax error tolerance to reduce required order
    devs = [0.01 0.005];  % relaxed error requirements
    
    % Design with kaiserord
    [n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fs);
    
    % Limit filter order to avoid excessive order
    if n > max_filter_order
        % If computed order too high, cap to a smaller order
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

% Method 2: Butterworth IIR filter function
function [filtered_signal, params] = apply_IIR_filter(signal, fs, cutoff_freq)
    % Design Butterworth IIR high-pass
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

% Method 3: Enhanced denoise + smoothing pre-pipeline (baseline removal + HP)
function [enhanced_signal, params] = apply_enhanced_pipeline(signal, fs)
    params = struct();

    % Dual-window median baseline removal (~0.2s and ~0.6s)
    win1 = max(3, 2*floor(0.20*fs)+1);
    win2 = max(3, 2*floor(0.60*fs)+1);
    if exist('medfilt1','file')
        baseline_rough = medfilt1(signal, win1, 'truncate');
        baseline = medfilt1(baseline_rough, win2, 'truncate');
    else
        % Fallback: moving median if medfilt1 not available
        baseline_rough = movmedian(signal, win1);
        baseline = movmedian(baseline_rough, win2);
    end
    signal_detrended = signal - baseline;

    % Mild high-pass to further suppress drift (no strict band-pass)
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