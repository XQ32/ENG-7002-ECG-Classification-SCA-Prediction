function view_shhs_ecg_v3()
% =================================================================================================
% File: view_shhs_ecg_v3.m
% Overview:
%   SHHS ECG PVC detection + record-level SCA risk quick assessment GUI (v3).
%   Single-axis scrollable window + multi-control workflow: Import EDF -> Preprocess -> Beat detection/classification -> PVC filter/navigation -> SCA risk prediction.
%
% Main responsibilities:
%   1. EDF import and ECG channel extraction (supports numeric/cell/string numeric) and basic metadata display (filename/size/duration).
%   2. Call ecgFilter (mode 2) to create filtered signals for visualization and detection (currently both use mode 2; reserved for future divergence).
%   3. Call detectAndClassifyHeartbeats to get R/QRS locations and beatInfo/stats.
%   4. Extract features using extractHeartbeatFeatures -> load beats classifier (results/trainedClassifier_latest.mat) -> predict PVC/Other.
%   5. Interactions:
%      - Time window slider browsing and second-level jumping
%      - Mouse wheel centered zoom / button zoom / reset
%      - PVC/Other/All filtering view
%      - Previous/next beat navigation + center on selection + highlight selected beat
%      - R point markers + labels (PVC red / Other gray)
%   6. SCA risk: use predicted PVC sequence -> aggregate record-level features on 5s post-PVC windows -> load results/sca_classifier_post_ectopic.mat -> output risk label.
%   7. Stats display: total beats, PVC/Other counts, selection index, risk result, and actual survival status (from CVD summary dataset).
%
% I/O:
%   External call: no input args; run directly to launch GUI.
%   User interactions: "Select EDF…", time slider, jump field, zoom in/out/reset, Find PVCs, filter dropdown, prev/next beat, Predict SCA Risk.
%   Data sources:
%     - EDF: shhs/polysomnography/edfs
%     - Beat classifier: results/trainedClassifier_latest.mat
%     - SCA risk model: results/sca_classifier_post_ectopic.mat
%     - Patient vital: shhs/datasets/shhs-cvd-summary-dataset-*.csv
%   Output: no files written (all shown in GUI); internal state kept in the state struct.
%
% Key implementation notes:
%   - Centralized state struct (raw/filtered ECG, beatInfo, predicted labels, filter masks, view window, models, risk result).
%   - Dual filtered channels reserved (ecgVisFiltered and ecgDetFiltered) for future differentiated params.
%   - Prediction pipeline: feature consistency check -> missing required features -> fill NaN with median -> call predictFcn.
%   - Custom 5s post-PVC aggregation: RR, HRT (TO/TS), QRS duration, R amplitude, HR-related, density, recovery, and HRV stats (mean/std/median/min/max).
%   - Actual vital lookup: parse trailing digits in EDF filename as nsrrid and match with CVD summary (vital 0/1).
%   - Center-on-selection + dynamic window label when zooming.
%
% Robustness & safeguards:
%   - All core steps (EDF read/filter/detect/classify/SCA predict) wrapped in try/catch with dialogs to avoid crashes.
%   - Channel search strategies (exact ECG / fuzzy ecg, ekg) and tolerant concatenation for cell/string numerics.
%   - Missing features / missing model files -> clear errordlg to stop the step.
%   - NaN feature imputation via sample median (fallback 0) to mitigate outliers.
%   - Window/index bounds protection: setViewStart clamping; interpY boundary clamp.
%   - Auto fallback to first filtered index if current selection becomes invalid.
%
% Performance notes:
%   - Local rendering per window: only samples within the current index range.
%   - Label count equals beats in window; can add thinning if needed (kept straightforward now).
%   - Mouse wheel scales the window length to avoid global recompute.
%
% Use cases / non-use cases:
%   Use: quick PVC QC for one record, PVC distribution exploration, instant SCA risk indication.
%   Not for: batch offline evaluation / training set processing (use script pipelines instead).
%
% Changelog:
%   2025-08-30 Added unified header comments (standardized); core logic unchanged.
%
% Usage:
%   >> view_shhs_ecg_v3
%   1) Click "Select EDF…" -> 2) "Find PVCs" -> 3) browse via filter/navigation -> 4) "Predict SCA Risk".
%
% Dependencies:
%   ecgFilter.m, detectAndClassifyHeartbeats.m, extractHeartbeatFeatures.m
%   results/trainedClassifier_latest.mat
%   results/sca_classifier_post_ectopic.mat
%   SHHS CVD summary dataset (vital field)
% =================================================================================================
    % 共享状态
    state = struct();
    state.fs = 125;                % default SHHS1 sampling rate
    state.windowSec = 30;          % browsing window length (seconds)
    state.viewStartSec = 0;        % current window start (sec)
    state.ecg = [];
    state.ecgFiltered = [];
    state.ecgVisFiltered = [];
    state.ecgDetFiltered = [];
    state.time = [];
    state.N = 0;
    state.edfPath = '';
    state.ecgVarName = '';
    state.fileSizeBytes = NaN;
    state.beatInfo = [];
    state.stats = struct();
    state.rGlobalAll = [];
    state.predLabels = {};
    state.pvcMask = [];
    state.otherMask = [];
    state.filteredBeatIndices = [];
    state.selectedBeatGlobalIdx = 1;
    state.filterMode = 'All'; % options: 'All'|'PVC'|'Other'
    state.clsModel = [];      % beat classifier model
    state.scaModel = [];      % SCA risk model
    state.scaResultText = '';
    state.patientVital = NaN; % 0=Dead,1=Alive

    % Build UI
    ui = buildUI();
    updateUIState('init');

    % Internal: build GUI
    function ui = buildUI()
        screenSize = get(0, 'ScreenSize');
        figW = min(1200, screenSize(3) * 0.9);
        figH = min(720, screenSize(4) * 0.85);
        ui.fig = figure('Name','SHHS ECG PVC & SCA Risk Analysis', 'NumberTitle','off', ...
            'MenuBar','none', 'ToolBar','none', 'Color','w', 'Units','pixels', ...
            'Position',[100 100 figW figH]);

        % Layout: main axes at top, controls at bottom
        ctrlH = 230;
        ui.ax = axes('Parent', ui.fig, 'Units','pixels', 'Position', [60 ctrlH+30 figW-100 figH-ctrlH-80]);
        set(ui.ax,'Color','w');
        grid(ui.ax,'on'); hold(ui.ax,'on');
        xlabel(ui.ax,'Time (s)'); ylabel(ui.ax,'ECG (a.u.)');

        % Controls panel
        panel = uipanel('Parent', ui.fig, 'Units','pixels', 'Position',[10 10 figW-20 ctrlH-20], 'BorderType','none', 'BackgroundColor','w');

        % Row 1: import & basic info
        y1 = ctrlH - 50;
        ui.btnLoad = uicontrol(panel,'Style','pushbutton','String','Select EDF…','Units','pixels', ...
            'Position',[10 y1 120 28],'Callback',@onLoadEDF,'BackgroundColor','w');
        ui.textFile = uicontrol(panel,'Style','text','String','File: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y1 460 22],'BackgroundColor','w');
        ui.textSize = uicontrol(panel,'Style','text','String','Size: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[610 y1 160 22],'BackgroundColor','w');
        ui.textFs   = uicontrol(panel,'Style','text','String','Sampling rate: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[780 y1 120 22],'BackgroundColor','w');
        ui.textDur  = uicontrol(panel,'Style','text','String','Duration: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[910 y1 200 22],'BackgroundColor','w');

        % Row 2: browsing & positioning
        y2 = ctrlH - 90;
        ui.textWindow = uicontrol(panel,'Style','text','String','Window: 0-30 s','HorizontalAlignment','left', ...
            'Units','pixels','Position',[10 y2 150 22],'BackgroundColor','w');
        ui.sliderTime = uicontrol(panel,'Style','slider','Units','pixels','Position',[100 y2 620 22], ...
            'Min',0,'Max',1,'Value',0,'Callback',@onSliderTime,'BackgroundColor','w');
        ui.textJump = uicontrol(panel,'Style','text','String','Go to (s):','HorizontalAlignment','left', ...
            'Units','pixels','Position',[740 y2 70 22],'BackgroundColor','w');
        ui.editJump = uicontrol(panel,'Style','edit','String','0','Units','pixels','Position',[810 y2-2 80 26],'BackgroundColor','w');
        ui.btnJump = uicontrol(panel,'Style','pushbutton','String','Go','Units','pixels', ...
            'Position',[900 y2-2 60 26],'Callback',@onJump,'BackgroundColor','w');
        ui.btnZoomIn = uicontrol(panel,'Style','pushbutton','String','Zoom In','Units','pixels', ...
            'Position',[970 y2-2 60 26],'Callback',@onZoomIn,'BackgroundColor','w');
        ui.btnZoomOut = uicontrol(panel,'Style','pushbutton','String','Zoom Out','Units','pixels', ...
            'Position',[1040 y2-2 60 26],'Callback',@onZoomOut,'BackgroundColor','w');
        ui.btnResetView = uicontrol(panel,'Style','pushbutton','String','Reset','Units','pixels', ...
            'Position',[1100 y2-2 60 26],'Callback',@onZoomReset,'BackgroundColor','w');

        % Row 3: PVC finding & stats
        y3 = ctrlH - 130;
        ui.btnFindPVC = uicontrol(panel,'Style','pushbutton','String','Find PVCs','Units','pixels', ...
            'Position',[10 y3 120 28],'Callback',@onFindPVC,'BackgroundColor','w');
        ui.textBeats = uicontrol(panel,'Style','text','String','Beats: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y3 200 22],'BackgroundColor','w');
        ui.textPVC   = uicontrol(panel,'Style','text','String','PVC: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[350 y3 160 22],'BackgroundColor','w');
        ui.textOther = uicontrol(panel,'Style','text','String','Other: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[520 y3 160 22],'BackgroundColor','w');

        % Row 4: navigation & filter
        y4 = ctrlH - 170;
        ui.btnPrev = uicontrol(panel,'Style','pushbutton','String','←','Units','pixels', ...
            'Position',[10 y4 60 36],'Callback',@(~,~) stepSelection(-1),'BackgroundColor','w','FontSize',12,'FontWeight','bold');
        ui.btnNext = uicontrol(panel,'Style','pushbutton','String','→','Units','pixels', ...
            'Position',[75 y4 60 36],'Callback',@(~,~) stepSelection(1),'BackgroundColor','w','FontSize',12,'FontWeight','bold');
        ui.textSel = uicontrol(panel,'Style','text','String','Selected: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[145 y4 220 26],'BackgroundColor','w','FontSize',10);
        ui.textFilter = uicontrol(panel,'Style','text','String','Filter:','HorizontalAlignment','left', ...
            'Units','pixels','Position',[370 y4 40 26],'BackgroundColor','w');
        ui.popupFilter = uicontrol(panel,'Style','popupmenu','String',{'All','PVC','Other'}, ...
            'Units','pixels','Position',[415 y4-2 120 28],'Callback',@onFilterChanged,'BackgroundColor','w');

        % Row 5: SCA risk
        y5 = ctrlH - 210;
        ui.btnSCA = uicontrol(panel,'Style','pushbutton','String','Predict SCA Risk','Units','pixels', ...
            'Position',[10 y5 120 28],'Callback',@onPredictSCA,'BackgroundColor','w');
        ui.textSCA = uicontrol(panel,'Style','text','String','SCA Risk: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y5 960 22],'BackgroundColor','w');

        % Pre-create graphics handles
        ui.hECG = plot(ui.ax, NaN, NaN, 'b-','LineWidth',1); % ECG trace
        ui.hR_PVC = plot(ui.ax, NaN, NaN, 'rv','MarkerFaceColor','r','MarkerSize',6,'LineStyle','none'); % PVC R
        ui.hR_Other = plot(ui.ax, NaN, NaN, 'k^','MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',5,'LineStyle','none'); % Other R
        ui.hSel = plot(ui.ax, NaN, NaN, 'go','MarkerSize',10,'LineWidth',2,'LineStyle','none'); % selected highlight
        ui.textGroup = gobjects(0); % text labels

        % Enable built-in zoom/pan and refresh
        ui.zoomObj = zoom(ui.fig); ui.zoomObj.Motion = 'both'; ui.zoomObj.Enable = 'on';
        ui.panObj = pan(ui.fig); ui.panObj.Enable = 'on';
        ui.zoomObj.ActionPostCallback = @onZoomPan;
        ui.panObj.ActionPostCallback = @onZoomPan;
        set(ui.fig,'WindowScrollWheelFcn', @onScrollZoom);
    end

    % Mouse wheel zoom (scale x-axis around cursor position)
    function onScrollZoom(~,evd)
        if isempty(state.ecg), return; end
        try
            cp = get(ui.ax,'CurrentPoint'); xCenter = cp(1,1);
        catch
            xl = get(ui.ax,'XLim'); xCenter = mean(xl);
        end
        xl = get(ui.ax,'XLim');
        xLen = max(1/state.fs, xl(2)-xl(1));
        if evd.VerticalScrollCount < 0
            scale = 0.8; % zoom in
        else
            scale = 1.25; % zoom out
        end
        newLen = max(1/state.fs, xLen*scale);
        left = xCenter - (xCenter - xl(1)) * (newLen/xLen);
        right = left + newLen;
        durSec = state.N / max(state.fs,1);
        if left < 0, left = 0; right = left + newLen; end
        if right > durSec, right = durSec; left = max(0, right - newLen); end
        set(ui.ax,'XLim',[left right]);
        set(ui.textWindow,'String',sprintf('Window: %.1f-%.1f s', left, right));
    end

    % Zoom in
    function onZoomIn(~,~)
        fakeEvt.VerticalScrollCount = -1; onScrollZoom([],fakeEvt);
    end
    % Zoom out
    function onZoomOut(~,~)
        fakeEvt.VerticalScrollCount = 1; onScrollZoom([],fakeEvt);
    end
    % Reset view
    function onZoomReset(~,~)
        setViewStart(0);
    end

    % Initialize or update UI state
    function updateUIState(~)
        % Basic control availability
        if isempty(state.ecg)
            set(ui.btnFindPVC,'Enable','off');
            set(ui.btnSCA,'Enable','off');
            set(ui.btnPrev,'Enable','off');
            set(ui.btnNext,'Enable','off');
            set(ui.popupFilter,'Enable','off');
            set(ui.sliderTime,'Enable','off');
            set(ui.btnJump,'Enable','off');
            set(ui.editJump,'Enable','off');
        else
            set(ui.btnFindPVC,'Enable','on');
            set(ui.sliderTime,'Enable','on');
            set(ui.btnJump,'Enable','on');
            set(ui.editJump,'Enable','on');
        end
        if ~isempty(state.rGlobalAll)
            set(ui.btnPrev,'Enable','on');
            set(ui.btnNext,'Enable','on');
            set(ui.popupFilter,'Enable','on');
            set(ui.btnSCA,'Enable','on');
        end

        % Text info
        if isempty(state.edfPath)
            set(ui.textFile,'String','File: -');
            set(ui.textSize,'String','Size: -');
            set(ui.textFs,'String','Sampling rate: -');
            set(ui.textDur,'String','Duration: -');
        else
            [~,nm,ext] = fileparts(state.edfPath);
            set(ui.textFile,'String',sprintf('File: %s%s', nm, ext));
            if isfinite(state.fileSizeBytes)
                set(ui.textSize,'String',sprintf('Size: %.2f MB', state.fileSizeBytes/1024/1024));
            else
                set(ui.textSize,'String','Size: -');
            end
            set(ui.textFs,'String',sprintf('Sampling rate: %d Hz', state.fs));
            durSec = state.N / max(state.fs,1);
            set(ui.textDur,'String',sprintf('Duration: %.1f s', durSec));
        end

        % Statistics
        if isempty(state.rGlobalAll)
            set(ui.textBeats,'String','Beats: -');
            set(ui.textPVC,'String','PVC: -');
            set(ui.textOther,'String','Other: -');
            set(ui.textSel,'String','Selected: -');
        else
            numBeats = numel(state.rGlobalAll);
            numPVC = sum(state.pvcMask);
            numOther = sum(state.otherMask);
            set(ui.textBeats,'String',sprintf('Beats: %d', numBeats));
            set(ui.textPVC,'String',sprintf('PVC: %d', numPVC));
            set(ui.textOther,'String',sprintf('Other: %d', numOther));
            % selection info
            sel = state.selectedBeatGlobalIdx;
            if ~isempty(sel) && sel>=1 && sel<=numBeats
                set(ui.textSel,'String',sprintf('Selected: %d/%d', findCurrentFilteredPos(), numel(state.filteredBeatIndices)));
            else
                set(ui.textSel,'String','Selected: -');
            end
        end

        % SCA result
        if isempty(state.scaResultText)
            set(ui.textSCA,'String','SCA Risk: -');
        else
            set(ui.textSCA,'String',state.scaResultText);
        end

        % Re-render
        renderECGWindow(true);
    end

    % Load EDF
    function onLoadEDF(~,~)
        try
            startDir = fullfile(pwd,'shhs','polysomnography','edfs');
            if ~isfolder(startDir), startDir = pwd; end
            [fn,fp] = uigetfile({'*.edf','EDF File (*.edf)'}, 'Select EDF file', startDir);
            if isequal(fn,0), return; end
            edfPath = fullfile(fp, fn);
            state.edfPath = edfPath;
            d = dir(edfPath); if ~isempty(d), state.fileSizeBytes = d.bytes; else, state.fileSizeBytes = NaN; end

            % Read EDF
            try
                TT = edfread(edfPath);
            catch ME
                errordlg(['edfread failed: ' ME.message], 'Read Error');
                return;
            end

            % Find ECG channel
            varNames = TT.Properties.VariableNames;
            ecgIdx = find(strcmp(varNames, 'ECG'), 1);
            if isempty(ecgIdx)
                vlow = lower(varNames);
                ecgIdx = find(contains(vlow, 'ecg') | contains(vlow, 'ekg'), 1, 'first');
            end
            if isempty(ecgIdx)
                errordlg('ECG channel not found.','Channel Error');
                return;
            end
            state.ecgVarName = varNames{ecgIdx};
            ecgCol = TT.(state.ecgVarName);
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
                errordlg(['Unsupported ECG channel type: ' class(ecgCol)], 'Channel Error');
                return;
            end
            state.ecg = double(ecg);
            state.N = numel(state.ecg);
            state.time = (0:state.N-1)'/state.fs;

            % Dual filtering: visualization (2), detection (2 here); SHHS1 has 60Hz notch
            try
                [state.ecgVisFiltered, ~] = ecgFilter(state.ecg, state.fs, 2, 0);
            catch ME
                warndlg(['Visualization filtering failed (using raw): ' ME.message], 'Filter Warning');
                state.ecgVisFiltered = state.ecg;
            end
            try
                [state.ecgDetFiltered, ~] = ecgFilter(state.ecg, state.fs, 2, 0);
            catch ME
                warndlg(['Detection filtering failed (fallback to visualization signal): ' ME.message], 'Filter Warning');
                state.ecgDetFiltered = state.ecgVisFiltered;
            end

            % Reset detection/results
            state.beatInfo = [];
            state.stats = struct();
            state.rGlobalAll = [];
            state.predLabels = {};
            state.pvcMask = [];
            state.otherMask = [];
            state.filteredBeatIndices = [];
            state.selectedBeatGlobalIdx = 1;
            state.filterMode = 'All';
            set(ui.popupFilter,'Value',1);
            state.scaResultText = '';

            % Reset view
            setViewStart(0);
            % Configure slider range
            durSec = state.N / max(state.fs,1);
            maxSlider = max(0, durSec - state.windowSec);
            set(ui.sliderTime,'Min',0,'Max',maxSlider,'Value',0);
            set(ui.textWindow,'String',sprintf('Window: %.1f-%.1f s', state.viewStartSec, state.viewStartSec+state.windowSec));

            updateUIState('loaded');
        catch ME
            errordlg(['Import failed: ' ME.message], 'Error');
        end
    end

    % Time slider
    function onSliderTime(~,~)
        val = get(ui.sliderTime,'Value');
        setViewStart(val);
    end

    % Quick jump
    function onJump(~,~)
        if isempty(state.ecg), return; end
        txt = get(ui.editJump,'String');
        t = str2double(txt);
        if ~isfinite(t), return; end
        durSec = state.N / max(state.fs,1);
        t = max(0, min(durSec - state.windowSec, t));
        set(ui.sliderTime,'Value',t);
        setViewStart(t);
    end

    % Set window start (sec) and refresh
    function setViewStart(t0)
        if isempty(state.ecg)
            state.viewStartSec = 0;
        else
            durSec = state.N / max(state.fs,1);
            t0 = max(0, min(max(0,durSec - state.windowSec), t0));
            state.viewStartSec = t0;
        end
        set(ui.textWindow,'String',sprintf('Window: %.1f-%.1f s', state.viewStartSec, state.viewStartSec+state.windowSec));
        renderECGWindow(false);
    end

    % Find PVCs
    function onFindPVC(~,~)
        if isempty(state.ecgDetFiltered)
            errordlg('Please import an EDF first.','Notice');
            return;
        end
        drawnow;
        try
            ATRTIMED = [];
            ANNOTD = {};
            [~, beatInfo, stats] = detectAndClassifyHeartbeats(state.ecgDetFiltered, ATRTIMED, ANNOTD, state.fs);
            if isfield(stats,'mteo_failed') && stats.mteo_failed
                errordlg('MTEO (Q/S) detection failed.','Detection Failed');
                return;
            end
            if isempty(beatInfo)
                errordlg('No valid beats found.','Detection Failed');
                return;
            end
            state.beatInfo = beatInfo;
            state.stats = stats;

            % Feature extraction
            [featureTable, ~] = extractHeartbeatFeatures(beatInfo, state.fs);
            if isempty(featureTable) || height(featureTable) == 0
                errordlg('Feature extraction returned empty.','Error');
                return;
            end

            % Load model
            modelFile = fullfile(pwd,'results','trainedClassifier_latest.mat');
            if ~exist(modelFile,'file')
                errordlg(['Model not found: ' modelFile],'Missing Model');
                return;
            end
            M = load(modelFile);
            if isfield(M,'trainedModelPackage')
                trainedClassifier = M.trainedModelPackage.trainedClassifier;
            else
                trainedClassifier = M.trainedClassifier;
            end
            state.clsModel = trainedClassifier;

            % Select variables and impute
            if isfield(trainedClassifier,'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
                reqVars = trainedClassifier.RequiredVariables;
            else
                reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
            end
            miss = setdiff(reqVars, featureTable.Properties.VariableNames);
            if ~isempty(miss)
                errordlg(['Feature table missing model required variables: ' strjoin(miss, ', ')],'Missing Variables');
                return;
            end
            X = featureTable(:, reqVars);
            for vv = 1:numel(reqVars)
                col = X.(reqVars{vv});
                if isnumeric(col)
                    nanMask = isnan(col);
                    if any(nanMask)
                        medv = median(col(~nanMask));
                        if isempty(medv) || isnan(medv), medv = 0; end
                        col(nanMask) = medv;
                        X.(reqVars{vv}) = col;
                    end
                end
            end

            % Predict
            try
                [predLabels, ~] = trainedClassifier.predictFcn(X);
            catch
                predLabels = trainedClassifier.predictFcn(X);
            end
            if ~iscell(predLabels), predLabels = cellstr(predLabels); end
            state.predLabels = predLabels;

            % Global R indices
            rGlobalAll = zeros(numel(beatInfo),1);
            for k = 1:numel(beatInfo)
                rGlobalAll(k) = double(beatInfo(k).segmentStartIndex) + double(beatInfo(k).rIndex) - 1;
            end
            state.rGlobalAll = rGlobalAll;

            % Masks
            state.pvcMask = strcmp(predLabels,'PVC');
            state.otherMask = strcmp(predLabels,'Other');

            % Filtering and selection
            rebuildFilteredIndices();
            if ~isempty(state.filteredBeatIndices)
                state.selectedBeatGlobalIdx = state.filteredBeatIndices(1);
                centerOnSelectedBeat();
            end

            updateUIState('predicted_pvc');
        catch ME
            errordlg(['Find PVCs failed: ' ME.message], 'Error');
        end
    end

    % Predict SCA risk (record-level aggregated features)
    function onPredictSCA(~,~)
        if isempty(state.rGlobalAll) || ~any(state.pvcMask)
            warndlg('No PVC beats available; cannot predict SCA risk.','Notice');
            return;
        end
        try
            % Load SCA model
            modelFile = fullfile(pwd,'results','sca_classifier_post_ectopic.mat');
            if ~exist(modelFile,'file')
                errordlg(['SCA model not found: ' modelFile],'Missing Model');
                return;
            end
            M = load(modelFile);
            scaModel = M.trainedClassifier;
            predictorNames = M.predictorNames;
            state.scaModel = scaModel;

            % Compute record-level aggregated features (aligned with training)
            pvcRidx = state.rGlobalAll(state.pvcMask);
            Trec = computePostPVCFeatureTable(pvcRidx, 5.0); % single-row aggregation table
            if isempty(Trec) || height(Trec)~=1
                errordlg('Did not obtain record-level aggregated features.','Error');
                return;
            end

            % Select model predictors and impute
            miss = setdiff(predictorNames, Trec.Properties.VariableNames);
            if ~isempty(miss)
                errordlg(['Missing required variables for SCA model: ' strjoin(miss, ', ')],'Missing Variables');
                return;
            end
            Xin = Trec(:, predictorNames);
            for i = 1:numel(predictorNames)
                col = Xin.(predictorNames{i});
                if ~isnumeric(col), continue; end
                nanMask = isnan(col);
                if any(nanMask)
                    medv = median(col(~nanMask));
                    if isempty(medv) || isnan(medv), medv = 0; end
                    col(nanMask) = medv;
                    Xin.(predictorNames{i}) = col;
                end
            end

            % Predict death risk probability (scores for class 0)
            try
                [~, rawScores] = predict(scaModel.ClassificationEnsemble, Xin);
            catch
                [~, rawScores] = predict(scaModel.ClassificationEnsemble, Xin{:,:});
            end
            try
                cls = scaModel.ClassificationEnsemble.ClassNames;
            catch
                cls = [0;1];
            end
            col0 = find(cls==0, 1, 'first'); if isempty(col0), col0=1; end
            scoreDead = rawScores(:, col0);

            thr = 0.5;
            if isfield(M,'modelOptions') && isfield(M.modelOptions,'enableThresholdMoving') && M.modelOptions.enableThresholdMoving
                if isfield(M.modelOptions,'pvcThreshold') && isfinite(M.modelOptions.pvcThreshold)
                    thr = M.modelOptions.pvcThreshold;
                end
            elseif isfield(scaModel,'Options') && isfield(scaModel.Options,'enableThresholdMoving') && scaModel.Options.enableThresholdMoving
                if isfield(scaModel.Options,'pvcThreshold') && isfinite(scaModel.Options.pvcThreshold)
                    thr = scaModel.Options.pvcThreshold;
                end
            end

            hasRisk = (scoreDead >= thr);

            % Lookup patient's actual vital
            [vitalVal, vitalFound] = getVitalForCurrentRecord();
            state.patientVital = vitalVal;
            vitalStr = 'Unknown';
            if vitalFound
                if isequal(vitalVal,1)
                    vitalStr = 'Alive(1)';
                elseif isequal(vitalVal,0)
                    vitalStr = 'Dead(0)';
                else
                    vitalStr = sprintf('%.0f', vitalVal);
                end
            end

            % pvcCount = Trec.record_pvc_count;
            % pvcPerHour = Trec.PVCs_per_hour;
            % state.scaResultText = sprintf('SCA Risk: %s | PVC count=%d, PVC density=%.2f /h | Death risk=%.3f, Thr=%.2f | Actual vital=%s', ...
            %     ternary(hasRisk,'At Risk','No Risk'), pvcCount, pvcPerHour, scoreDead, thr, vitalStr);
            state.scaResultText = sprintf('SCA Risk: %s | Actual vital=%s', ...
            ternary(hasRisk,'At Risk','No Risk'), vitalStr);

            updateUIState('sca_predicted');
        catch ME
            errordlg(['SCA risk prediction failed: ' ME.message],'Error');
        end
    end

    % Render current window
    function renderECGWindow(preserveXLim)
        if nargin<1, preserveXLim=false; end
        if isempty(state.ecg)
            set(ui.hECG,'XData',NaN,'YData',NaN);
            set(ui.hR_PVC,'XData',NaN,'YData',NaN);
            set(ui.hR_Other,'XData',NaN,'YData',NaN);
            set(ui.hSel,'XData',NaN,'YData',NaN);
            deleteTexts();
            return;
        end
        fs = state.fs;
        N = state.N;
        if preserveXLim
            xl = get(ui.ax,'XLim');
            t0 = max(0, xl(1));
            w = max(1/fs, xl(2) - xl(1));
        else
            t0 = state.viewStartSec;
            w = state.windowSec;
        end
        iStart = max(1, floor(t0*fs)+1);
        iEnd = min(N, iStart + round(w*fs) - 1);
        t = (iStart:iEnd)'/fs;
        if ~isempty(state.ecgVisFiltered)
            y = state.ecgVisFiltered(iStart:iEnd);
        else
            y = state.ecg(iStart:iEnd);
        end
        set(ui.hECG,'XData',t,'YData',y);
        if ~preserveXLim
            xlim(ui.ax,[t0 t0+w]);
        end
        % R markers
        deleteTexts();
        if ~isempty(state.rGlobalAll)
            inView = state.rGlobalAll>=iStart & state.rGlobalAll<=iEnd;
            rIdxView = state.rGlobalAll(inView);
            tR = (rIdxView - 1)/fs;
            lblView = state.predLabels(inView);
            isPVC = strcmp(lblView,'PVC');
            isOther = strcmp(lblView,'Other');
            % Positions from detection indices, amplitude from visualization series
            set(ui.hR_PVC,'XData',tR(isPVC),'YData',interpY(iStart, y, rIdxView(isPVC)));
            set(ui.hR_Other,'XData',tR(isOther),'YData',interpY(iStart, y, rIdxView(isOther)));
            % Text labels
            ui.textGroup = gobjects(sum(inView),1);
            rYRange = max(y) - min(y);
            for k = 1:numel(rIdxView)
                tx = tR(k); ty = interpY(iStart, y, rIdxView(k)) + 0.05*rYRange;
                col = [0.2 0.2 0.2];
                if strcmp(lblView{k},'PVC'), col = [0.85 0.1 0.1]; end
                ui.textGroup(k) = text(ui.ax, tx, ty, lblView{k}, 'Color', col, 'FontSize',8, 'HorizontalAlignment','center');
            end
            % Selection highlight
            sel = state.selectedBeatGlobalIdx;
            if ~isempty(sel) && sel>=1 && sel<=numel(state.rGlobalAll)
                selIdx = state.rGlobalAll(sel);
                if selIdx>=iStart && selIdx<=iEnd
                    tSel = (selIdx-1)/fs;
                    ySel = interpY(iStart, y, selIdx);
                    set(ui.hSel,'XData',tSel,'YData',ySel);
                else
                    set(ui.hSel,'XData',NaN,'YData',NaN);
                end
            else
                set(ui.hSel,'XData',NaN,'YData',NaN);
            end
        else
            set(ui.hR_PVC,'XData',NaN,'YData',NaN);
            set(ui.hR_Other,'XData',NaN,'YData',NaN);
            set(ui.hSel,'XData',NaN,'YData',NaN);
        end
        drawnow;
    end

    % Refresh after zoom/pan
    function onZoomPan(~,~)
        if isempty(state.ecg), return; end
        renderECGWindow(true);
        try
            xl = get(ui.ax,'XLim');
            set(ui.textWindow,'String',sprintf('Window: %.1f-%.1f s', xl(1), xl(2)));
        catch
        end
    end

    % Delete text labels
    function deleteTexts()
        if ~isempty(ui.textGroup)
            mask = isgraphics(ui.textGroup);
            if any(mask)
                delete(ui.textGroup(mask));
            end
        end
        ui.textGroup = gobjects(0);
    end

    % Interpolate y at index (clamped to bounds)
    function yv = interpY(iStart, yseg, globalIdx)
        % globalIdx is into the original signal; convert to segment index
        ii = max(1, min(numel(yseg), globalIdx - iStart + 1));
        yv = yseg(ii);
    end

    % Rebuild filtered indices
    function rebuildFilteredIndices()
        if isempty(state.rGlobalAll)
            state.filteredBeatIndices = [];
            return;
        end
        switch state.filterMode
            case 'PVC'
                mask = state.pvcMask;
            case 'Other'
                mask = state.otherMask;
            otherwise
                mask = true(size(state.rGlobalAll));
        end
        state.filteredBeatIndices = find(mask(:));
        if isempty(state.filteredBeatIndices)
            state.filteredBeatIndices = [];
        end
    end

    % Position of current selection among filtered indices
    function pos = findCurrentFilteredPos()
        pos = NaN;
        if isempty(state.filteredBeatIndices), return; end
        pos = find(state.filteredBeatIndices == state.selectedBeatGlobalIdx, 1, 'first');
        if isempty(pos)
            % if not found, find the nearest
            dif = abs(state.filteredBeatIndices - state.selectedBeatGlobalIdx);
            [~,mi] = min(dif);
            pos = mi;
        end
    end

    % Move selection left/right
    function stepSelection(step)
        if isempty(state.filteredBeatIndices), return; end
        curPos = findCurrentFilteredPos();
        if ~isfinite(curPos), curPos = 1; end
        newPos = max(1, min(numel(state.filteredBeatIndices), curPos + step));
        state.selectedBeatGlobalIdx = state.filteredBeatIndices(newPos);
        centerOnSelectedBeat();
        updateUIState('step_sel');
    end

    % Filter changed
    function onFilterChanged(~,~)
        items = get(ui.popupFilter,'String');
        val = get(ui.popupFilter,'Value');
        state.filterMode = items{val};
        rebuildFilteredIndices();
        if ~isempty(state.filteredBeatIndices)
            state.selectedBeatGlobalIdx = state.filteredBeatIndices(1);
            centerOnSelectedBeat();
        end
        updateUIState('filter_changed');
    end

    % Center view on selected beat
    function centerOnSelectedBeat()
        if isempty(state.rGlobalAll) || isempty(state.selectedBeatGlobalIdx), return; end
        fs = state.fs; w = state.windowSec;
        selIdx = state.rGlobalAll(state.selectedBeatGlobalIdx);
        tSel = (selIdx-1)/fs;
        durSec = state.N / max(fs,1);
        t0 = max(0, min(max(0,durSec - w), tSel - w/2));
        set(ui.sliderTime,'Value',t0);
        setViewStart(t0);
    end

    % Compute record-level aggregated features (aligned with training script: variable observation window, recovery, clusters/density, etc.)
    function T = computePostPVCFeatureTable(pvcIndices, postWindowSec)
        if isempty(pvcIndices)
            T = table(); return;
        end
        fs = state.fs; N = state.N;
        rGlobalAll = state.rGlobalAll(:);
        beatInfo = state.beatInfo;

        % Config (consistent with training script)
        baselineSec = 40; baselineMinBeats = 10; maxObsSec = 20; consecBeats = 5;
        rrTolFrac = 0.10; rrTolSigma = 2.0; tsLowThr = 0.0; hrTolBpm = 5;
        joinSec = 10; winSec30 = 30; winSec60 = 60;

        % Sort PVC global samples and map to beat indices
        pvcIndicesSorted = sort(pvcIndices(:));
        idxPVC_all_sorted = round(interp1(rGlobalAll, 1:numel(rGlobalAll), pvcIndicesSorted, 'nearest', 'extrap'));
        idxPVC_all_sorted = max(1, min(numel(rGlobalAll), idxPVC_all_sorted));
        numPVC = numel(pvcIndicesSorted);

        % RR vectors
        numBeatsRec = numel(rGlobalAll);
        rr_between = nan(numBeatsRec,1);
        if numBeatsRec >= 2
            rr_between(2:end) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
        end
        rr_pre_vec = nan(numBeatsRec,1); rr_post1_vec = nan(numBeatsRec,1); rr_post2_vec = nan(numBeatsRec,1);
        if numBeatsRec >= 2
            rr_pre_vec(2:end) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
            rr_post1_vec(1:end-1) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
        end
        if numBeatsRec >= 3
            rr_post2_vec(1:end-2) = (rGlobalAll(3:end) - rGlobalAll(1:end-2)) / fs;
        end

        % QRS duration and R amplitude
        qrs_dur_vec = nan(numBeatsRec,1); r_amp_vec = nan(numBeatsRec,1);
        for ii = 1:numBeatsRec
            b = beatInfo(ii);
            on = NaN; off = NaN;
            if isfield(b,'qrsOnIndex'), on = double(b.qrsOnIndex); end
            if isfield(b,'qrsOffIndex'), off = double(b.qrsOffIndex); end
            if isfinite(on) && isfinite(off) && off > on
                qrs_dur_vec(ii) = (off - on) / fs;
            end
            if isfield(b,'segment') && isfield(b,'rIndex') && ~isempty(b.segment) && isfinite(b.rIndex) && b.rIndex>=1 && b.rIndex<=numel(b.segment)
                r_amp_vec(ii) = double(b.segment(b.rIndex));
            end
        end

        % PVC flags and SQI
        isPVCBeat = false(numBeatsRec,1);
        if ~isempty(idxPVC_all_sorted)
            isPVCBeat(unique(idxPVC_all_sorted)) = true;
        end
        if isfield(beatInfo, 'sqiIsGood')
            try
                sqiVec = arrayfun(@(b) (isfield(b,'sqiIsGood') && ~isempty(b.sqiIsGood) && logical(b.sqiIsGood)), beatInfo);
                sqiVec = logical(sqiVec(:));
            catch
                sqiVec = true(numBeatsRec,1);
            end
        else
            sqiVec = true(numBeatsRec,1);
        end

        % Per-PVC measures: recovery, AUC, HRT, TS abnormality
        nextPVCsamples = [pvcIndicesSorted(2:end); inf];
        recovered_flags = false(numPVC,1);
        recovery_time_sec = nan(numPVC,1);
        auc_rr_dev = nan(numPVC,1);
        hrt_to_vals = nan(numPVC,1);
        hrt_ts_abnormal = nan(numPVC,1);

        for kk = 1:numPVC
            pvcSample = pvcIndicesSorted(kk);
            idxPVC = idxPVC_all_sorted(kk);
            nextPVC = nextPVCsamples(kk);
            obsEnd = min([double(N), double(pvcSample) + round(maxObsSec*fs), double(nextPVC)-1]);
            if ~isfinite(obsEnd) || obsEnd <= pvcSample
                continue;
            end

            % Personalized baseline RR
            baseA = max(1, double(pvcSample) - round(baselineSec*fs));
            baseB = max(1, double(pvcSample) - 1);
            j_in_baseline = find(rGlobalAll >= baseA & rGlobalAll <= baseB);
            j_in_baseline = j_in_baseline(:);
            j_in_baseline = j_in_baseline(j_in_baseline >= 2);
            if ~isempty(j_in_baseline)
                mask_ok = ~isPVCBeat(j_in_baseline) & ~isPVCBeat(j_in_baseline-1) & sqiVec(j_in_baseline) & sqiVec(j_in_baseline-1);
                j_in_baseline = j_in_baseline(mask_ok);
            end
            rr_base_vec = rr_between(j_in_baseline);
            rr_base_vec = rr_base_vec(isfinite(rr_base_vec));
            if numel(rr_base_vec) < max(3, baselineMinBeats)
                j_all = (2:numBeatsRec).';
                mask_all = ~isPVCBeat(j_all) & ~isPVCBeat(j_all-1) & sqiVec(j_all) & sqiVec(j_all-1);
                rr_all = rr_between(j_all(mask_all));
                rr_all = rr_all(isfinite(rr_all));
                if isempty(rr_all)
                    muRR = median(rr_between(isfinite(rr_between)));
                    sigRR = std(rr_between(isfinite(rr_between)));
                else
                    muRR = mean(rr_all);
                    sigRR = std(rr_all);
                end
            else
                muRR = mean(rr_base_vec);
                sigRR = std(rr_base_vec);
            end
            if ~isfinite(muRR) || muRR<=0
                muRR = median(rr_between(isfinite(rr_between)));
            end
            if ~isfinite(sigRR)
                sigRR = 0.0;
            end

            % HRT (TO/TS) and abnormality
            rr_pre = NaN; rr_post1 = NaN; rr_post2 = NaN;
            if idxPVC >= 2 && isfinite(rr_between(idxPVC)), rr_pre = rr_between(idxPVC); end
            if idxPVC+1 <= numBeatsRec && isfinite(rr_between(idxPVC+1)), rr_post1 = rr_between(idxPVC+1); end
            if idxPVC+2 <= numBeatsRec && isfinite(rr_between(idxPVC+2)), rr_post2 = rr_between(idxPVC+2); end
            if isfinite(rr_pre) && isfinite(rr_post1)
                hrt_to_vals(kk) = (rr_post1 - rr_pre) / max(rr_pre, eps);
            else
                hrt_to_vals(kk) = NaN;
            end
            % TS: maximum average slope of any 5 consecutive non-PVC good-SQI RR within 5-15 beats after PVC
            j_after = (idxPVC+1):min(numBeatsRec, idxPVC+20); j_after = j_after(:);
            mask_rr_ok = true(size(j_after));
            mask_rr_ok = mask_rr_ok & (rGlobalAll(j_after) <= obsEnd) & (rGlobalAll(j_after-1) >= pvcSample);
            mask_rr_ok = mask_rr_ok & ~isPVCBeat(j_after) & ~isPVCBeat(j_after-1) & sqiVec(j_after) & sqiVec(j_after-1);
            j_after = j_after(mask_rr_ok);
            ts_val = NaN;
            if numel(j_after) >= 7
                rr_seq = rr_between(j_after);
                maxSlope = -inf;
                for uu = 1:(numel(rr_seq)-4)
                    slope = (rr_seq(uu+4) - rr_seq(uu)) / 4.0;
                    if slope > maxSlope, maxSlope = slope; end
                end
                ts_val = maxSlope;
            end
            if isfinite(ts_val)
                hrt_ts_abnormal(kk) = double(ts_val <= tsLowThr);
            else
                hrt_ts_abnormal(kk) = NaN;
            end

            % Recovery (RR-based) and AUC(|RR-μRR|)
            j = idxPVC + 1; consec = 0; recFlag = false; recJ = NaN; aucVal = 0.0; lastSample = pvcSample;
            while j <= numBeatsRec
                if rGlobalAll(j) > obsEnd, break; end
                if j >= 2 && ~isPVCBeat(j) && ~isPVCBeat(j-1) && sqiVec(j) && sqiVec(j-1) && isfinite(rr_between(j))
                    rrj = rr_between(j);
                    tolAbs = max(rrTolFrac * muRR, rrTolSigma * sigRR);
                    if abs(rrj - muRR) <= tolAbs
                        consec = consec + 1;
                        if consec >= consecBeats
                            recFlag = true; recJ = j; 
                            % 仍然累计AUC到此点
                        end
                    else
                        consec = 0;
                    end
                    dt = (rGlobalAll(j) - lastSample) / baselineSec; % time normalization
                    aucVal = aucVal + abs(rrj - muRR) * max(dt, eps);
                    lastSample = rGlobalAll(j);
                    if recFlag, break; end
                end
                j = j + 1;
            end
            recovered_flags(kk) = recFlag;
            if recFlag && isfinite(recJ)
                recovery_time_sec(kk) = (double(rGlobalAll(recJ)) - double(pvcSample)) / fs;
            else
                recovery_time_sec(kk) = (double(obsEnd) - double(pvcSample)) / fs;
            end
            auc_rr_dev(kk) = aucVal;
        end

        % Assemble window-level feature vectors (consistent with training)
        RR_Pre = rr_pre_vec(idxPVC_all_sorted);
        RR_Post1 = rr_post1_vec(idxPVC_all_sorted);
        RR_Post2 = rr_post2_vec(idxPVC_all_sorted);
        HRT_TO = (RR_Post1 - RR_Pre) ./ (RR_Pre + eps);
        HRT_TS = (RR_Post2 - RR_Pre) / 2;
        QRS_Dur_PVC = qrs_dur_vec(idxPVC_all_sorted);
        R_Amp_PVC = r_amp_vec(idxPVC_all_sorted);
        % Keep legacy names for compatibility; placeholders retained
        Beats_in_5s = nan(numPVC,1); HR_5s = nan(numPVC,1);
        HR_Pre = 60 ./ RR_Pre; HR_Post1 = 60 ./ RR_Post1;

        % Derived vectors
        HR_Accel = HR_Post1 - HR_Pre;
        CompRatio = RR_Post1 ./ (RR_Pre + eps);
        pvcIntervals = diff(pvcIndicesSorted) / fs;

        % Clusters and density
        pvcTimes = pvcIndicesSorted / fs;
        clusterStarts = []; clusterEnds = []; clusterSizes = [];
        if numPVC >= 1
            curStart = pvcTimes(1); curEnd = pvcTimes(1); curSize = 1;
            for kk = 2:numPVC
                if (pvcTimes(kk) - curEnd) < joinSec
                    curEnd = pvcTimes(kk); curSize = curSize + 1;
                else
                    clusterStarts(end+1) = curStart; %#ok<AGROW>
                    clusterEnds(end+1) = curEnd; %#ok<AGROW>
                    clusterSizes(end+1) = curSize; %#ok<AGROW>
                    curStart = pvcTimes(kk); curEnd = pvcTimes(kk); curSize = 1;
                end
            end
            clusterStarts(end+1) = curStart; clusterEnds(end+1) = curEnd; clusterSizes(end+1) = curSize;
        end
        clusterDurations = max(0, clusterEnds - clusterStarts);
        pvcDensity30_max_per_min = max_density(pvcTimes, winSec30);
        pvcDensity60_max_per_min = max_density(pvcTimes, winSec60);
        longest_noPVC_sec = NaN;
        if numPVC >= 2
            longest_noPVC_sec = max(diff(pvcTimes));
        elseif numPVC == 1
            longest_noPVC_sec = max(pvcTimes(1), (N/fs) - pvcTimes(1));
        end

        % Record-level rates
        recordDurationHr = (N / fs) / 3600;
        pvcPerHour = numPVC / max(recordDurationHr, eps);

        % Build aggregated structure (field names aligned with training)
        S = struct();
        S.record_pvc_count = numPVC;
        S.PVCs_per_hour = pvcPerHour;

        S = agg_add_stats(S, RR_Pre,    'RR_Pre');
        S = agg_add_stats(S, RR_Post1,  'RR_Post1');
        S = agg_add_stats(S, RR_Post2,  'RR_Post2');
        S = agg_add_stats(S, HRT_TO,    'HRT_TO');
        S = agg_add_stats(S, HRT_TS,    'HRT_TS');
        S = agg_add_stats(S, QRS_Dur_PVC, 'QRS_Dur_PVC');
        S = agg_add_stats(S, R_Amp_PVC, 'R_Amp_PVC');
        S = agg_add_stats(S, Beats_in_5s, 'Beats_in_5s');
        S = agg_add_stats(S, HR_Pre,    'HR_Pre');
        S = agg_add_stats(S, HR_Post1,  'HR_Post1');
        S = agg_add_stats(S, HR_5s,     'HR_5s');

        % Additional aggregates
        mask_to = isfinite(HRT_TO);
        if any(mask_to), S.HRT_TO_neg_frac = mean(double(HRT_TO(mask_to) < 0)); else, S.HRT_TO_neg_frac = NaN; end
        mask_qrs = isfinite(QRS_Dur_PVC);
        if any(mask_qrs), S.QRS_Prolonged_frac = mean(double(QRS_Dur_PVC(mask_qrs) > 0.12)); else, S.QRS_Prolonged_frac = NaN; end
        S = agg_add_stats(S, HR_Accel,   'HR_Accel');
        S = agg_add_stats(S, CompRatio,  'CompRatio');
        S.RR_Pre_CV   = agg_cv(RR_Pre);
        S.RR_Post1_CV = agg_cv(RR_Post1);
        S.RR_Pre_RMSSD   = agg_rmssd(RR_Pre);
        S.RR_Post1_RMSSD = agg_rmssd(RR_Post1);
        S = agg_add_stats(S, recovery_time_sec, 'recovery_time_sec');
        S = agg_add_stats(S, auc_rr_dev, 'AUC_RR_Deviation');
        S = agg_add_stats(S, pvcIntervals, 'PVC_Interval');

        % Clusters & high density
        S.cluster_count = numel(clusterSizes);
        S = agg_add_stats(S, clusterDurations, 'cluster_duration_sec');
        S = agg_add_stats(S, clusterSizes,     'cluster_size');
        S.PVC_density30_max_per_min = pvcDensity30_max_per_min;
        S.PVC_density60_max_per_min = pvcDensity60_max_per_min;
        S.longest_noPVC_sec = longest_noPVC_sec;

        % Recovery ratios
        S.recovered_frac = mean(double(recovered_flags));
        S.recovery_failure_ratio = mean(double(~recovered_flags));

        % HRT abnormal fractions
        if any(isfinite(hrt_to_vals))
            S.HRT_TO_abnormal_frac = mean(double(hrt_to_vals(isfinite(hrt_to_vals)) >= 0));
        else
            S.HRT_TO_abnormal_frac = NaN;
        end
        if any(isfinite(hrt_ts_abnormal))
            S.HRT_TS_abnormal_frac = mean(double(hrt_ts_abnormal(isfinite(hrt_ts_abnormal)) > 0));
        else
            S.HRT_TS_abnormal_frac = NaN;
        end

        % Recovery failure fraction inside clusters
        if ~isempty(clusterSizes)
            rec_fail_in_clusters = [];
            pvcTimesVec = pvcTimes; %#ok<NASGU>
            for c = 1:numel(clusterSizes)
                inC = pvcTimes >= clusterStarts(c) & pvcTimes <= clusterEnds(c);
                if any(inC)
                    rec_fail_in_clusters(end+1) = mean(double(~recovered_flags(inC))); %#ok<AGROW>
                end
            end
            if ~isempty(rec_fail_in_clusters)
                S.recovery_failure_in_clusters_frac_mean = mean(rec_fail_in_clusters);
                S.recovery_failure_in_clusters_frac_max  = max(rec_fail_in_clusters);
            else
                S.recovery_failure_in_clusters_frac_mean = NaN;
                S.recovery_failure_in_clusters_frac_max  = NaN;
            end
        else
            S.recovery_failure_in_clusters_frac_mean = NaN;
            S.recovery_failure_in_clusters_frac_max  = NaN;
        end

        T = struct2table(S);

        % Local density function
        function mx = max_density(ts, winSec)
            if isempty(ts)
                mx = 0; return;
            end
            i = 1; j = 0; mx = 0; nloc = numel(ts);
            for i = 1:nloc
                while j < nloc && (ts(j+1) - ts(i)) <= winSec
                    j = j + 1;
                end
                if (j - i + 1) > mx
                    mx = j - i + 1;
                end
            end
            mx = mx / (winSec/60);
        end
    end

    % Read patient's actual vital (vital: 0=Dead,1=Alive)
    function [vitalVal, found] = getVitalForCurrentRecord()
        vitalVal = NaN; found = false;
        if isempty(state.edfPath), return; end
        [~, base, ~] = fileparts(state.edfPath);
        tok = regexp(base, '([0-9]+)$', 'tokens', 'once');
        if isempty(tok)
            tok = regexp(base, '([0-9]{6,})', 'tokens', 'once');
        end
        if isempty(tok), return; end
        idStr = tok{1};
        idNum = str2double(idStr);
        % Candidate paths
        p1 = fullfile(pwd,'shhs','datasets','shhs-cvd-summary-dataset-0.21.0.csv');
        cand = p1;
        if ~exist(cand,'file')
            d1 = dir(fullfile(pwd,'shhs','datasets','**','shhs-cvd-summary-dataset-*.csv'));
            if ~isempty(d1)
                % pick the most recently modified
                [~,mi] = max([d1.datenum]); %#ok<DATNM>
                cand = fullfile(d1(mi).folder, d1(mi).name);
            else
                cand = '';
            end
        end
        if isempty(cand) || ~exist(cand,'file'), return; end
        try
            T = readtable(cand);
        catch
            return;
        end
        vn = lower(T.Properties.VariableNames);
        iId = find(strcmp(vn,'nsrrid'),1); if isempty(iId), iId = find(contains(vn,'nsrrid'),1); end
        iVital = find(strcmp(vn,'vital'),1); if isempty(iVital), iVital = find(contains(vn,'vital'),1); end
        if isempty(iId) || isempty(iVital), return; end
        idCol = T{:,iId};
        if isnumeric(idCol)
            ids = double(idCol);
        else
            ids = str2double(string(idCol));
        end
        vitalCol = T{:,iVital};
        vtmp = string(vitalCol);
        % match
        idx = find(ids == idNum, 1, 'first');
        if isempty(idx), return; end
        % parse vital
        valStr = lower(strtrim(vtmp(idx)));
        if any(strcmp(valStr, {'1','alive'}))
            vitalVal = 1; found = true; return;
        elseif any(strcmp(valStr, {'0','dead'}))
            vitalVal = 0; found = true; return;
        else
            vNum = str2double(valStr);
            if isfinite(vNum)
                vitalVal = vNum; found = true; return;
            end
        end
    end

    % Ternary helper
    function out = ternary(cond, a, b)
        if cond, out = a; else, out = b; end
    end

    % ====== Record-level aggregation helpers ======
    function S = agg_add_stats(S, v, baseName)
        if nargin < 3 || isempty(baseName)
            baseName = 'feat';
        end
        if isempty(v)
            vals = [NaN, NaN, NaN, NaN, NaN];
        else
            v = v(:);
            v = v(isfinite(v));
            if isempty(v)
                vals = [NaN, NaN, NaN, NaN, NaN];
            else
                m = mean(v);
                s = std(v);
                med = median(v);
                mn = min(v);
                mx = max(v);
                vals = [m, s, med, mn, mx];
            end
        end
        S.([baseName '_mean']) = vals(1);
        S.([baseName '_std']) = vals(2);
        S.([baseName '_median']) = vals(3);
        S.([baseName '_min']) = vals(4);
        S.([baseName '_max']) = vals(5);
    end

    function cv = agg_cv(v)
        if isempty(v)
            cv = NaN; return; end
        v = v(:);
        v = v(isfinite(v));
        if numel(v) < 2
            cv = NaN; return; end
        mu = mean(v);
        sigma = std(v);
        den = abs(mu);
        if den < 1e-12
            cv = NaN; return; end
        cv = sigma / den;
    end

    function rmssd = agg_rmssd(v)
        if isempty(v) || numel(v) < 2
            rmssd = NaN; return; end
        v = v(:);
        v = v(isfinite(v));
        if numel(v) < 2
            rmssd = NaN; return; end
        d = diff(v);
        rmssd = sqrt(mean(d.^2));
    end
end
