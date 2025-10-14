function view_shhs_ecg_v3()
% =================================================================================================
% File: view_shhs_ecg_v3.m
% Overview:
%   SHHS ECG PVC detection + record-level SCA risk rapid assessment integrated GUI (version 3).
%   Single-axis scrolling window + multi-control operations: Import EDF -> Preprocess -> Beat detection/classification -> PVC filtering/navigation -> SCA risk prediction.
%
% Main responsibilities:
%   1. EDF import and ECG channel extraction (compatible with numeric / cell / numeric string) and basic metadata display (file name/size/duration).
%   2. Call ecgFilter (mode 2) to generate filtered signals for visualization and detection (currently both visualization/detection use mode 2, with room for differentiation).
%   3. Call detectAndClassifyHeartbeats to perform R/QRS localization and obtain basic beatInfo/statistics.
%   4. Based on extractHeartbeatFeatures features -> load beats classification model (results/trainedClassifier_latest.mat) -> predict PVC/Other.
%   5. Interactive features:
%      - Time window slider browsing and second-level jump
%      - Mouse wheel centered zoom / button zoom / reset
%      - Filter views: PVC/Other/All
%      - Previous/next beat navigation + center positioning + selected beat highlight
%      - R point markers + label text (PVC red / Other gray)
%   6. SCA risk: Use predicted PVC sequence -> aggregate record-level features from 5-second windows -> load results/sca_classifier_post_ectopic.mat classifier -> output risk label.
%   7. Statistics display: total beats, PVC/Other counts, selected position, risk result and actual vital status (datasets CVD summary).
%
% I/O:
%   External call: No input parameters; run directly to start the GUI.
%   User interactions: "Choose EDF…", time slider, jump input, zoom in/out/reset buttons, find PVC, filter dropdown, previous/next beat, predict SCA risk.
%   Data sources:
%     - EDF: shhs/polysomnography/edfs
%     - Beat model: results/trainedClassifier_latest.mat
%     - SCA risk model: results/sca_classifier_post_ectopic.mat
%     - Patient vital info: shhs/datasets/shhs-cvd-summary-dataset-*.csv
%   Output: No file write (all displayed in GUI); internal state stored in the state struct.
%
% Key implementation points:
%   - Centralized state struct (raw/filtered ECG, beatInfo, predicted labels, filter mask, view window, models, risk results).
%   - Dual-filter channel design reserved (ecgVisFiltered and ecgDetFiltered) to allow parameter differentiation later.
%   - Prediction flow: feature consistency check -> reject missing features -> fill NaN with sample median (fallback 0 when empty) -> call predictFcn.
%   - Custom 5-second post-PVC aggregation: RR, HRT (TO/TS), QRS duration, R amplitude, HR related, density, recovery and HRV statistics (mean/std/median/min/max).
%   - Actual vital info lookup: parse trailing digits of EDF filename as nsrrid and match with CVD summary table (vital 0/1).
%   - View center positioning + dynamic window label updates after zooming.
%
% Robustness and safeguards:
%   - All core steps (EDF read/filter/detect/classify/SCA predict) wrapped in try/catch with dialogs/warnings to avoid crashes.
%   - Multi-strategy channel search (exact ECG / fuzzy ecg, ekg) and tolerant concatenation for cell/numeric string.
%   - Missing features / missing model files -> explicit errordlg to stop corresponding flow.
%   - NaN feature imputation uses sample median (fallback 0) to reduce outlier impact.
%   - Window boundary/index bounds protection: setViewStart with clamping; interpY shrinks to boundaries.
%   - If selection becomes invalid after filtering, automatically fallback to the first filtered index.
%
% Performance notes:
%   - Local rendering per window: only take current window index range (based on start second and length).
%   - Number of labels equals number of beats in the window; if further optimization needed, down-sampling can be added (kept intuitive for now).
%   - Wheel zoom scales window length proportionally to avoid global recomputation.
%
% Applicable/Not applicable:
%   Applicable: Quick PVC QC for a single record, PVC distribution exploration, instant view of potential SCA risk indication.
%   Not applicable: Large-scale offline evaluation / batch training set processing (use scripting pipelines).
%
% Change log:
%   2025-08-30 Unified header comments added (standardized Chinese structure; core logic unchanged).
%
% Usage:
%   >> view_shhs_ecg_v3
%   1) Click "Choose EDF…" -> 2) "Find PVC" -> 3) Browse via filter/navigation -> 4) "Predict SCA Risk".
%
% Dependencies:
%   ecgFilter.m, detectAndClassifyHeartbeats.m, extractHeartbeatFeatures.m
%   results/trainedClassifier_latest.mat
%   results/sca_classifier_post_ectopic.mat
%   SHHS CVD summary dataset (vital field)
% =================================================================================================
    % Shared state
    state = struct();
    state.fs = 125;                % Default SHHS1 sampling rate
    state.windowSec = 30;          % View window length (seconds)
    state.viewStartSec = 0;        % Current window start second
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
    state.filterMode = 'All'; % Options: 'All'|'PVC'|'Other'
    state.clsModel = [];      % Beat classification model
    state.scaModel = [];      % SCA risk model
    state.scaResultText = '';
    state.patientVital = NaN; % 0=Dead,1=Alive

    % Build UI
    ui = buildUI();
    updateUIState('init');

    % Internal: build interface
    function ui = buildUI()
        screenSize = get(0, 'ScreenSize');
        figW = min(1200, screenSize(3) * 0.9);
        figH = min(720, screenSize(4) * 0.85);
        ui.fig = figure('Name','SHHS ECG PVC & SCA Risk Analysis', 'NumberTitle','off', ...
            'MenuBar','none', 'ToolBar','none', 'Color','w', 'Units','pixels', ...
            'Position',[100 100 figW figH]);

        % Layout: large axis on top, controls at bottom
        ctrlH = 230;
        ui.ax = axes('Parent', ui.fig, 'Units','pixels', 'Position', [60 ctrlH+30 figW-100 figH-ctrlH-80]);
        set(ui.ax,'Color','w');
        grid(ui.ax,'on'); hold(ui.ax,'on');
        xlabel(ui.ax,'Time (s)'); ylabel(ui.ax,'ECG (a.u.)');

        % Control area (using pixels layout)
        panel = uipanel('Parent', ui.fig, 'Units','pixels', 'Position',[10 10 figW-20 ctrlH-20], 'BorderType','none', 'BackgroundColor','w');

        % Row 1: import and basic info
        y1 = ctrlH - 50;
        ui.btnLoad = uicontrol(panel,'Style','pushbutton','String','Choose EDF…','Units','pixels', ...
            'Position',[10 y1 120 28],'Callback',@onLoadEDF,'BackgroundColor','w');
        ui.textFile = uicontrol(panel,'Style','text','String','File: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y1 460 22],'BackgroundColor','w');
        ui.textSize = uicontrol(panel,'Style','text','String','Size: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[610 y1 160 22],'BackgroundColor','w');
        ui.textFs   = uicontrol(panel,'Style','text','String','Sample rate: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[780 y1 120 22],'BackgroundColor','w');
        ui.textDur  = uicontrol(panel,'Style','text','String','Duration: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[910 y1 200 22],'BackgroundColor','w');

        % Row 2: browsing and positioning
        y2 = ctrlH - 90;
        ui.textWindow = uicontrol(panel,'Style','text','String','Window: 0-30 s','HorizontalAlignment','left', ...
            'Units','pixels','Position',[10 y2 150 22],'BackgroundColor','w');
        ui.sliderTime = uicontrol(panel,'Style','slider','Units','pixels','Position',[100 y2 620 22], ...
            'Min',0,'Max',1,'Value',0,'Callback',@onSliderTime,'BackgroundColor','w');
        ui.textJump = uicontrol(panel,'Style','text','String','Jump to (s):','HorizontalAlignment','left', ...
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

        % Row 3: PVC finding and stats
        y3 = ctrlH - 130;
        ui.btnFindPVC = uicontrol(panel,'Style','pushbutton','String','Find PVC','Units','pixels', ...
            'Position',[10 y3 120 28],'Callback',@onFindPVC,'BackgroundColor','w');
        ui.textBeats = uicontrol(panel,'Style','text','String','Beats: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y3 200 22],'BackgroundColor','w');
        ui.textPVC   = uicontrol(panel,'Style','text','String','PVC: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[350 y3 160 22],'BackgroundColor','w');
        ui.textOther = uicontrol(panel,'Style','text','String','Other: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[520 y3 160 22],'BackgroundColor','w');

        % Row 4: navigation and filter
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

        % Pre-created plotting handles
        ui.hECG = plot(ui.ax, NaN, NaN, 'b-','LineWidth',1); % ECG main line
        ui.hR_PVC = plot(ui.ax, NaN, NaN, 'rv','MarkerFaceColor','r','MarkerSize',6,'LineStyle','none'); % PVC R
        ui.hR_Other = plot(ui.ax, NaN, NaN, 'k^','MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',5,'LineStyle','none'); % Other R
        ui.hSel = plot(ui.ax, NaN, NaN, 'go','MarkerSize',10,'LineWidth',2,'LineStyle','none'); % Selected highlight
        ui.textGroup = gobjects(0); % Text labels

        % Enable built-in zoom/pan and link refresh
        ui.zoomObj = zoom(ui.fig); ui.zoomObj.Motion = 'both'; ui.zoomObj.Enable = 'on';
        ui.panObj = pan(ui.fig); ui.panObj.Enable = 'on';
        ui.zoomObj.ActionPostCallback = @onZoomPan;
        ui.panObj.ActionPostCallback = @onZoomPan;
        set(ui.fig,'WindowScrollWheelFcn', @onScrollZoom);
    end

    % Mouse wheel zoom (zoom x-axis around cursor)
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
        % Basic button enable state
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

        % Text information
        if isempty(state.edfPath)
            set(ui.textFile,'String','File: -');
            set(ui.textSize,'String','Size: -');
            set(ui.textFs,'String','Sample rate: -');
            set(ui.textDur,'String','Duration: -');
        else
            [~,nm,ext] = fileparts(state.edfPath);
            set(ui.textFile,'String',sprintf('File: %s%s', nm, ext));
            if isfinite(state.fileSizeBytes)
                set(ui.textSize,'String',sprintf('Size: %.2f MB', state.fileSizeBytes/1024/1024));
            else
                set(ui.textSize,'String','Size: -');
            end
            set(ui.textFs,'String',sprintf('Sample rate: %d Hz', state.fs));
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
            % Selected info
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

        % Plot refresh
        renderECGWindow(true);
    end

    % Load EDF
    function onLoadEDF(~,~)
        try
            startDir = fullfile(pwd,'shhs','polysomnography','edfs');
            if ~isfolder(startDir), startDir = pwd; end
            [fn,fp] = uigetfile({'*.edf','EDF Files (*.edf)'}, 'Select EDF file', startDir);
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

            % Dual filtering: visualization (2), detection (2); SHHS1 already has 60Hz notch
            try
                [state.ecgVisFiltered, ~] = ecgFilter(state.ecg, state.fs, 2, 0);
            catch ME
                warndlg(['Visualization filtering failed (using raw signal): ' ME.message], 'Filter Warning');
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

    % Set window start second and refresh
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

    % Find PVC
    function onFindPVC(~,~)
        if isempty(state.ecgDetFiltered)
            errordlg('Please load an EDF first.','Info');
            return;
        end
        drawnow;
        try
            ATRTIMED = [];
            ANNOTD = {};
            [~, beatInfo, stats] = detectAndClassifyHeartbeats(state.ecgDetFiltered, ATRTIMED, ANNOTD, state.fs);
            if isfield(stats,'mteo_failed') && stats.mteo_failed
                errordlg('MTEO(Q/S) detection failed.','Detection Failed');
                return;
            end
            if isempty(beatInfo)
                errordlg('No valid heartbeats detected.','Detection Failed');
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
                errordlg(['Model not found: ' modelFile],'Model Missing');
                return;
            end
            M = load(modelFile);
            if isfield(M,'trainedModelPackage')
                trainedClassifier = M.trainedModelPackage.trainedClassifier;
            else
                trainedClassifier = M.trainedClassifier;
            end
            state.clsModel = trainedClassifier;

            % Select variables and impute missing
            if isfield(trainedClassifier,'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
                reqVars = trainedClassifier.RequiredVariables;
            else
                reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
            end
            miss = setdiff(reqVars, featureTable.Properties.VariableNames);
            if ~isempty(miss)
                errordlg(['Features missing required variables for model: ' strjoin(miss, ', ')],'Missing Variables');
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
            errordlg(['Find PVC pipeline failed: ' ME.message], 'Error');
        end
    end

    % Predict SCA risk (record-level aggregated features)
    function onPredictSCA(~,~)
        if isempty(state.rGlobalAll) || ~any(state.pvcMask)
            warndlg('No PVC beats available; cannot predict SCA risk.','Info');
            return;
        end
        try
            % Load SCA model (v3)
            modelFile = fullfile(pwd,'results','sca_classifier_v3.mat');
            if ~exist(modelFile,'file')
                errordlg(['SCA model not found: ' modelFile],'Model Missing');
                return;
            end
            M = load(modelFile);
            scaModel = M.trainedClassifier;
            predictorNames = M.predictorNames;
            state.scaModel = scaModel;

            % Compute record-level aggregated features (aligned with training)
            pvcRidx = state.rGlobalAll(state.pvcMask);
            Trec = computePostPVCFeatureTable(pvcRidx, 5.0); % single-row aggregated table
            if isempty(Trec) || height(Trec)~=1
                errordlg('Failed to obtain record-level aggregated features.','Error');
                return;
            end

            % Inject demographics/follow-up info (age_s1, gender, race, censdate), align to latest feature names
            try
                demo = getSummaryFieldsForCurrentRecord();
                fns = fieldnames(demo);
                for ii = 1:numel(fns)
                    fn = fns{ii};
                    if ~ismember(fn, Trec.Properties.VariableNames)
                        Trec.(fn) = demo.(fn);
                    else
                        if ~isfinite(Trec.(fn)) && isfinite(demo.(fn))
                            Trec.(fn) = demo.(fn);
                        end
                    end
                end
            catch
            end

            % Align with model required features: rename/derive/fill
            try
                Trec = alignFeaturesForSCA(Trec, predictorNames);
            catch
            end

            % Select required variables and fill missing
            miss = setdiff(predictorNames, Trec.Properties.VariableNames);
            if ~isempty(miss)
                % Best-effort completion: fields that cannot be computed set to 0
                for jj = 1:numel(miss)
                    try
                        Trec.(miss{jj}) = 0;
                    catch
                    end
                end
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

            % Predict death risk score (scores is P(y==0))
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

            % Query patient vital info
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
            % state.scaResultText = sprintf('SCA Risk: %s | PVC count=%d, PVC density=%.2f /h | Death risk=%.3f, Threshold=%.2f | Actual vital=%s', ...
            %     ternary(hasRisk,'At Risk','No Risk'), pvcCount, pvcPerHour, scoreDead, thr, vitalStr);
            state.scaResultText = sprintf('SCA Risk: %s | Actual Vital Status=%s', ...
            ternary(hasRisk,'At Risk','No Risk'), vitalStr);

            updateUIState('sca_predicted');
        catch ME
            errordlg(['SCA risk prediction failed: ' ME.message],'Error');
        end
    end

    % View rendering
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
            % Use R indices from detection filtered signal for location, but amplitudes from visualization filtered signal
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
            % Selected highlight
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

    % Interpolate y at selected point (avoid out-of-bounds)
    function yv = interpY(iStart, yseg, globalIdx)
        % globalIdx is original index, convert to segment index
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

    % Find current selection position within filtered sequence
    function pos = findCurrentFilteredPos()
        pos = NaN;
        if isempty(state.filteredBeatIndices), return; end
        pos = find(state.filteredBeatIndices == state.selectedBeatGlobalIdx, 1, 'first');
        if isempty(pos)
            % If not inside, find nearest
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

    % Compute record-level aggregated features (v3 version: recovery dynamics feature system, aligned with generate_sca_classifier_v3.m)
    function T = computePostPVCFeatureTable(pvcIndices, ~)
        if isempty(pvcIndices)
            T = table(); return;
        end
        fs = state.fs; N = state.N;
        rGlobalAll = state.rGlobalAll(:);
        beatInfo = state.beatInfo;

        % Config parameters (aligned with generate_sca_classifier_v3.m lines 75-86)
        baselineSec = 50;        % Baseline window before PVC (seconds)
        baselineMinBeats = 10;   % Minimum baseline valid sinus beats
        maxObsSec = 50;          % Max observation window after PVC (seconds)
        consecStableBeats = 10;  % Continuous stable beats required to consider recovered
        rrTolFrac = 0.08;        % RR deviation threshold fraction
        rrTolSigma = 1.5;        % RR deviation threshold in sigmas
        hrt_ts_low_thr = 0.0;    % HRT TS low threshold
        minRR = 0.30; maxRR = 2.50;  % Reasonable RR range

        % PVC global sample point sorting and mapping to beat index
        pvcIndicesSorted = sort(pvcIndices(:));
        idxPVC_all_sorted = round(interp1(rGlobalAll, 1:numel(rGlobalAll), pvcIndicesSorted, 'nearest', 'extrap'));
        idxPVC_all_sorted = max(1, min(numel(rGlobalAll), idxPVC_all_sorted));
        numPVC = numel(pvcIndicesSorted);

        % RR vector & base signal vectors
        numBeatsRec = numel(rGlobalAll);
        rr_between = nan(numBeatsRec,1);
        if numBeatsRec >= 2
            rr_between(2:end) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
        end
        
        % QRS and R amplitude vectors (for possible QTc, etc.)
        qrs_dur_vec = nan(numBeatsRec,1); 
        r_amp_vec = nan(numBeatsRec,1);
        t_amp_vec = nan(numBeatsRec,1);  % T-wave amplitude (if available)
        tGlobalIdx = nan(numBeatsRec,1); % T-wave index (if available)
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
            % T-wave info (optional)
            if isfield(b,'tAmplitude') && isfinite(b.tAmplitude)
                t_amp_vec(ii) = double(b.tAmplitude);
            end
            if isfield(b,'tIndex') && isfinite(b.tIndex)
                tGlobalIdx(ii) = double(b.segmentStartIndex) + double(b.tIndex) - 1;
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

        % Per-PVC new feature system (aligned with generate_sca_classifier_v3.m lines 297-503)
        % Initialize feature containers
        halflife_rr       = nan(numPVC,1);  % 50% deviation recovery half-life
        halflife_rr_30    = nan(numPVC,1);  % 30% deviation recovery half-life
        halflife_rr_10    = nan(numPVC,1);  % 10% deviation recovery half-life
        tau_rr            = nan(numPVC,1);  % Recovery time constant
        auc_rr            = nan(numPVC,1);  % Recovery curve AUC
        overshoot_rr      = nan(numPVC,1);  % Overshoot magnitude
        oscill_rr         = nan(numPVC,1);  % Oscillation index
        stalls_rr         = nan(numPVC,1);  % Number of recovery stalls
        notRecovered_rr   = false(numPVC,1);% Not recovered flag
        
        hrt_to            = nan(numPVC,1);  % HRT turbulence onset
        hrt_ts            = nan(numPVC,1);  % HRT turbulence slope
        coupling_ratio    = nan(numPVC,1);  % Coupling interval ratio
        hr_peak_accel_bpm = nan(numPVC,1);  % Peak heart rate acceleration
        late_var_rr       = nan(numPVC,1);  % Late RR variability
        pvc_to_next_nonPVC_sec = nan(numPVC,1); % PVC to next non-PVC interval
        
        % High-discrimination/complexity features
        early_mean_dev_rr = nan(numPVC,1);  % Early mean deviation
        late_mean_dev_rr  = nan(numPVC,1);  % Late mean deviation
        slope_rr_dev_0_5s = nan(numPVC,1);  % 0-5s slope
        slope_rr_dev_5_15s= nan(numPVC,1);  % 5-15s slope
        sampen_rr         = nan(numPVC,1);  % Sample entropy
        lzc_rr            = nan(numPVC,1);  % Lempel-Ziv complexity
        poincare_ratio_pp = nan(numPVC,1);  % Poincaré ratio
        dev_rr_q95        = nan(numPVC,1);  % 95th percentile deviation
        tamp_rr_corr      = nan(numPVC,1);  % T-wave - RR correlation
        
        % T-wave and QTc recovery (simplified, GUI may not have full data)
        halflife_tamp     = nan(numPVC,1);
        auc_tamp          = nan(numPVC,1);
        halflife_qtc      = nan(numPVC,1);
        auc_qtc           = nan(numPVC,1);
        qtc_over_frac     = nan(numPVC,1);
        qtc_over_mag      = nan(numPVC,1);

        % PVC loop
        for kk = 1:numPVC
            pvcSample = pvcIndicesSorted(kk);
            idxPVC = idxPVC_all_sorted(kk);
            
            nextPVCsample = inf;
            if kk < numPVC
                nextPVCsample = pvcIndicesSorted(kk+1);
            end
            obsEnd = min([double(N), double(pvcSample)+round(maxObsSec*fs), double(nextPVCsample)-1]);
            if ~(isfinite(obsEnd) && obsEnd > pvcSample)
                notRecovered_rr(kk) = true;
                continue;
            end

            % Individualized baseline (aligned with generate_sca_classifier_v3.m lines 348-353 local_baseline_stats)
            baseA = max(1, double(pvcSample) - round(baselineSec*fs));
            baseB = max(1, double(pvcSample) - 1);
            j_in = find(rGlobalAll >= baseA & rGlobalAll <= baseB);
            j_in = j_in(:); j_in = j_in(j_in>=2);
            baseRR_vec = [];
            if ~isempty(j_in)
                mask_ok = ~isPVCBeat(j_in) & ~isPVCBeat(j_in-1) & sqiVec(j_in) & sqiVec(j_in-1) & isfinite(rr_between(j_in));
                baseRR_vec = rr_between(j_in(mask_ok));
            end
            baseRR_vec = baseRR_vec(isfinite(baseRR_vec) & baseRR_vec>=minRR & baseRR_vec<=maxRR);
            
            if numel(baseRR_vec) < max(3, baselineMinBeats)
                % Fallback to entire-record sinus RR
                j_all = (2:numBeatsRec).';
                mask_all = ~isPVCBeat(j_all) & ~isPVCBeat(j_all-1) & sqiVec(j_all) & sqiVec(j_all-1) & isfinite(rr_between(j_all));
                rr_all = rr_between(j_all(mask_all)); rr_all = rr_all(isfinite(rr_all));
                if isempty(rr_all)
                    muRR = median(rr_between(isfinite(rr_between)));
                    sigRR= std(rr_between(isfinite(rr_between)));
                else
                    muRR = median(rr_all);
                    sigRR= std(rr_all);
                end
            else
                muRR = mean(baseRR_vec);
                sigRR= std(baseRR_vec);
            end
            if ~isfinite(muRR) || muRR<=0
                muRR = median(rr_between(isfinite(rr_between)));
            end
            if ~isfinite(sigRR), sigRR = 0; end

            % Pre/post RR for HRT (aligned with generate_sca_classifier_v3.m lines 356-376)
            rr_pre=nan; rr_post1=nan; rr_post2=nan;
            if idxPVC>=2,          rr_pre   = rr_between(idxPVC); end
            if idxPVC+1<=numBeatsRec, rr_post1 = rr_between(idxPVC+1); end
            if idxPVC+2<=numBeatsRec, rr_post2 = rr_between(idxPVC+2); end
            
            if isfinite(rr_pre) && isfinite(rr_post1) && rr_pre>0
                hrt_to(kk) = (rr_post1 - rr_pre) / rr_pre;
            end
            % TS: maximum 5-consecutive slope
            hrt_ts(kk) = local_approx_ts(idxPVC, pvcSample, obsEnd, rGlobalAll, rr_between, isPVCBeat, sqiVec, numBeatsRec);
            
            if isfinite(rr_pre), coupling_ratio(kk) = rr_pre / muRR; end
            try
                hr0 = 60/max(muRR,eps);
                hr1 = 60/max(rr_post1,eps);
                hr2 = 60/max(rr_post2,eps);
                hr_peak_accel_bpm(kk) = max([hr1-hr0, hr2-hr0]);
            catch
                hr_peak_accel_bpm(kk)=NaN;
            end
            
            % Post-PVC non-PVC + good SQI sequence
            j_after = (idxPVC+1):numBeatsRec;
            keep = false(size(j_after));
            for jj = 1:numel(j_after)
                j = j_after(jj);
                if rGlobalAll(j) > obsEnd, break; end
                if j>=2 && ~isPVCBeat(j) && ~isPVCBeat(j-1) && sqiVec(j) && sqiVec(j-1) && ...
                        isfinite(rr_between(j)) && rr_between(j)>=minRR && rr_between(j)<=maxRR
                    keep(jj)=true;
                end
            end
            j_after = j_after(keep);
            if isempty(j_after)
                notRecovered_rr(kk)=true;
                continue;
            end
            t_after = (double(rGlobalAll(j_after)) - double(pvcSample))/fs;
            rr_seq  = rr_between(j_after);
            dev_rr  = abs(rr_seq - muRR)./max(muRR,eps);
            
            % Half-lives (50%/30%/10%)
            halflife_rr(kk)    = local_halflife(dev_rr, t_after, 0.50, consecStableBeats);
            halflife_rr_30(kk) = local_halflife(dev_rr, t_after, 0.30, consecStableBeats);
            halflife_rr_10(kk) = local_halflife(dev_rr, t_after, 0.10, consecStableBeats);
            if ~isfinite(halflife_rr(kk)) && t_after(end) < maxObsSec
                notRecovered_rr(kk) = true;
            end
            
            tau_rr(kk)      = local_time_constant(dev_rr, t_after, 8);
            auc_rr(kk)      = local_auc(dev_rr, t_after, baselineSec);
            overshoot_rr(kk)= max(dev_rr)-dev_rr(1);
            oscill_rr(kk)   = local_oscillation_index(dev_rr);
            stalls_rr(kk)   = local_recovery_stalls(dev_rr, 0.20, consecStableBeats);
            
            late_mask = t_after >= 5;
            if any(late_mask)
                late_var_rr(kk) = var(dev_rr(late_mask), 'omitnan');
            end
            
            % Complexity/slope/entropy/LZ (aligned with generate_sca_classifier_v3.m lines 436-483)
            early_mask = (t_after>=0)&(t_after<=5);
            late_mask2= (t_after>5)&(t_after<=15);
            if any(early_mask)
                early_mean_dev_rr(kk) = mean(dev_rr(early_mask),'omitnan');
                slope_rr_dev_0_5s(kk) = local_linear_slope(t_after(early_mask), dev_rr(early_mask));
            end
            if any(late_mask2)
                late_mean_dev_rr(kk)  = mean(dev_rr(late_mask2),'omitnan');
                slope_rr_dev_5_15s(kk)= local_linear_slope(t_after(late_mask2), dev_rr(late_mask2));
            end
            try
                sampen_rr(kk) = local_sampen(dev_rr, 2, 0.2*std(dev_rr,'omitnan'));
            catch
                sampen_rr(kk)=NaN;
            end
            try
                d1 = diff(dev_rr(:));
                lzc_rr(kk) = local_lzc_binary(double(d1>=0));
            catch
                lzc_rr(kk)=NaN;
            end
            try
                [SD1_pre, SD2_pre] = local_poincare_sd(baseRR_vec);
                post_mask10 = t_after <= 10;
                if any(post_mask10)
                    [SD1_post, SD2_post] = local_poincare_sd(rr_seq(post_mask10));
                    r_pre  = SD1_pre/max(SD2_pre,eps);
                    r_post = SD1_post/max(SD2_post,eps);
                    if isfinite(r_pre) && r_pre>0 && isfinite(r_post)
                        poincare_ratio_pp(kk) = r_post / r_pre;
                    end
                end
            catch
                poincare_ratio_pp(kk)=NaN;
            end
            try
                dev_rr_q95(kk) = quantile(dev_rr,0.95);
            catch
                dev_rr_q95(kk)=NaN;
            end
            
            % PVC → next non-PVC interval
            try
                j_search = (idxPVC+1):numBeatsRec;
                next_nonPVC_j = NaN;
                for jj = j_search
                    if ~isPVCBeat(jj)
                        next_nonPVC_j = jj; break;
                    end
                end
                if isfinite(next_nonPVC_j)
                    pvc_to_next_nonPVC_sec(kk) = (double(rGlobalAll(next_nonPVC_j)) - double(pvcSample))/fs;
                    if ~(isfinite(pvc_to_next_nonPVC_sec(kk)) && pvc_to_next_nonPVC_sec(kk)>=0)
                        pvc_to_next_nonPVC_sec(kk)=NaN;
                    end
                end
            catch
                pvc_to_next_nonPVC_sec(kk)=NaN;
            end
        end  % end PVC loop

        % ========== Record-level aggregation (aligned with generate_sca_classifier_v3.m lines 505-751) ==========
        recDurationHr = (N/fs)/3600;
        pvcPerHour = numPVC / max(recDurationHr,eps);
        recovery_failure_ratio = mean(double(notRecovered_rr));
        
        agg_mean = @(v) mean(v(isfinite(v)),'omitnan');
        agg_med  = @(v) median(v(isfinite(v)));
        agg_max  = @(v) local_safe_max(v);
        agg_min  = @(v) local_safe_min(v);
        
        R = struct();
        % Record name
        try
            [~,nm,~] = fileparts(state.edfPath);
            R.record = string(nm);
        catch
            R.record = "";
        end
        R.patientVital               = NaN;  % Fill later from CSV
        R.PVCs_per_hour              = pvcPerHour;
        R.PVC_count                  = numPVC;
        R.recovery_failure_ratio     = recovery_failure_ratio;
        R.halflife_rr_median         = agg_med(halflife_rr);
        R.halflife_rr30_median       = agg_med(halflife_rr_30);
        R.halflife_rr10_median       = agg_med(halflife_rr_10);
        R.tau_rr_median              = agg_med(tau_rr);
        R.oscill_rr_median           = agg_med(oscill_rr);
        R.HRT_TO_abnormal_frac       = mean(double(isfinite(hrt_to) & (hrt_to>=0)),'omitnan');
        R.HRT_TS_low_frac            = mean(double(isfinite(hrt_ts) & (hrt_ts<=hrt_ts_low_thr)),'omitnan');
        R.halflife_tamp_median       = agg_med(halflife_tamp);
        R.halflife_qtc_median        = agg_med(halflife_qtc);
        R.auc_qtc_mean               = agg_mean(auc_qtc);
        
        % RMSSD baseline (simplified)
        try
            R.RR_Pre_RMSSD_mean = NaN;  % Simplified in GUI
        catch
            R.RR_Pre_RMSSD_mean = NaN;
        end
        
        R.coupling_ratio_mean        = agg_mean(coupling_ratio);
        R.hr_peak_accel_bpm_mean     = agg_mean(hr_peak_accel_bpm);
        
        medEarly = agg_med(early_mean_dev_rr);
        medLate  = agg_med(late_mean_dev_rr);
        R.rr_dev_late_minus_early    = medLate - medEarly;
        R.rr_dev_slope_0_5s_med      = agg_med(slope_rr_dev_0_5s);
        R.tamp_rr_corr_med           = agg_med(tamp_rr_corr);
        R.qtc_overshoot_frac_mean    = agg_mean(qtc_over_frac);
        
        try
            hl = halflife_rr(isfinite(halflife_rr) & halflife_rr>0);
            if ~isempty(hl)
                R.recovery_rate_hmean = mean(1./hl);
            else
                R.recovery_rate_hmean = NaN;
            end
        catch
            R.recovery_rate_hmean = NaN;
        end
        
        R.PVC_to_next_nonPVC_sec_mean = agg_mean(pvc_to_next_nonPVC_sec);
        R.PVC_to_next_nonPVC_sec_max  = agg_max(pvc_to_next_nonPVC_sec);
        R.PVC_to_next_nonPVC_sec_min  = agg_min(pvc_to_next_nonPVC_sec);
        
        % ========== Innovative features (lines 566-751) ==========
        % 1. Nonlinear transforms
        R.oscill_log            = log(max(R.oscill_rr_median, eps));
        R.pvc_freq_sqrt         = sqrt(max(R.PVCs_per_hour, 0));
        R.halflife10_sqrt       = sqrt(max(R.halflife_rr10_median, eps));
        
        % 2. Feature interactions
        R.oscill_x_pvc_freq     = R.oscill_rr_median * R.PVCs_per_hour;
        R.oscill_x_recovery     = R.oscill_rr_median * R.halflife_rr10_median;
        R.pvc_x_variability     = R.PVCs_per_hour * R.RR_Pre_RMSSD_mean;
        
        % 3. Ratios
        R.oscill_per_pvc        = R.oscill_rr_median / max(R.PVCs_per_hour, eps);
        R.recovery_efficiency   = 1.0 / max(R.halflife_rr10_median, eps);
        R.early_late_ratio      = R.halflife_rr10_median / max(R.halflife_rr30_median, eps);
        
        % 4. Variability/stability
        try
            oscill_valid = oscill_rr(isfinite(oscill_rr) & oscill_rr > 0);
            if numel(oscill_valid) >= 3
                R.oscill_cv = std(oscill_valid) / max(mean(oscill_valid), eps);
            else
                R.oscill_cv = NaN;
            end
            
            hl_valid = halflife_rr(isfinite(halflife_rr) & halflife_rr > 0);
            if numel(hl_valid) >= 3
                R.halflife_cv = std(hl_valid) / max(mean(hl_valid), eps);
            else
                R.halflife_cv = NaN;
            end
            
            hl10_valid = halflife_rr(isfinite(halflife_rr));
            R.fast_recovery_ratio = mean(double(hl10_valid <= 10));
        catch
            R.oscill_cv = NaN;
            R.halflife_cv = NaN;
            R.fast_recovery_ratio = NaN;
        end
        
        % 5. Composite risk score
        try
            oscill_norm = min(R.oscill_rr_median / 0.5, 1.0);
            pvc_norm    = min(R.PVCs_per_hour / 50, 1.0);
            hl10_norm   = min(R.halflife_rr10_median / 30, 1.0);
            R.composite_risk_score = 0.40*oscill_norm + 0.30*pvc_norm + 0.30*hl10_norm;
        catch
            R.composite_risk_score = NaN;
        end
        
        % 6. Temporal dynamics
        try
            if isfinite(R.halflife_rr10_median) && isfinite(R.halflife_rr30_median) && ...
               R.halflife_rr10_median > 0 && R.halflife_rr30_median > R.halflife_rr10_median
                speed_10 = 0.10 / R.halflife_rr10_median;
                speed_30 = 0.20 / (R.halflife_rr30_median - R.halflife_rr10_median);
                R.recovery_deceleration = speed_10 - speed_30;
            else
                R.recovery_deceleration = NaN;
            end
        catch
            R.recovery_deceleration = NaN;
        end
        
        % 7. HRT composite
        try
            R.HRT_abnormal_combined = R.HRT_TO_abnormal_frac + R.HRT_TS_low_frac;
            if R.HRT_TO_abnormal_frac > 0.5 && R.HRT_TS_low_frac > 0.5
                R.HRT_risk_category = 2;
            elseif R.HRT_TO_abnormal_frac > 0.3 || R.HRT_TS_low_frac > 0.3
                R.HRT_risk_category = 1;
            else
                R.HRT_risk_category = 0;
            end
        catch
            R.HRT_abnormal_combined = NaN;
            R.HRT_risk_category = NaN;
        end
        
        % 8. More log transforms
        try
            R.pvc_freq_log = log(max(R.PVCs_per_hour, 1));
            R.halflife30_log = log(max(R.halflife_rr30_median, eps));
        catch
            R.pvc_freq_log = NaN;
            R.halflife30_log = NaN;
        end
        
        % 9. CV combination
        try
            R.combined_cv = (R.halflife_cv + R.oscill_cv) / 2;
            R.recovery_oscill_cv_ratio = R.halflife_cv / max(R.oscill_cv, eps);
        catch
            R.combined_cv = NaN;
            R.recovery_oscill_cv_ratio = NaN;
        end
        
        % 10. More interactions
        try
            R.oscill_log_x_pvc = R.oscill_log * R.PVCs_per_hour;
            R.oscill_x_cv = R.oscill_rr_median * R.halflife_cv;
            R.early_late_x_cv = R.early_late_ratio * R.halflife_cv;
        catch
            R.oscill_log_x_pvc = NaN;
            R.oscill_x_cv = NaN;
            R.early_late_x_cv = NaN;
        end
        
        % 11. Clinical composite indices
        try
            instability_score = 0.4*R.oscill_cv + 0.4*R.halflife_cv + 0.2*min(R.recovery_failure_ratio, 1.0);
            R.instability_score = instability_score;
            
            pvc_burden = (R.PVCs_per_hour / 50) * (R.oscill_rr_median / 0.5);
            R.pvc_burden_index = min(pvc_burden, 2.0);
            
            if isfinite(R.fast_recovery_ratio) && isfinite(R.early_late_ratio) && R.early_late_ratio > 0
                R.recovery_capacity = R.fast_recovery_ratio / R.early_late_ratio;
            else
                R.recovery_capacity = NaN;
            end
        catch
            R.instability_score = NaN;
            R.pvc_burden_index = NaN;
            R.recovery_capacity = NaN;
        end
        
        % 12. Ratios from top features
        try
            if isfinite(R.oscill_log) && isfinite(R.pvc_freq_log)
                R.oscill_pvc_log_ratio = R.oscill_log / max(R.pvc_freq_log, eps);
            else
                R.oscill_pvc_log_ratio = NaN;
            end
            R.tau_x_cv = R.tau_rr_median * R.halflife_cv;
        catch
            R.oscill_pvc_log_ratio = NaN;
            R.tau_x_cv = NaN;
        end
        
        % 13. Extreme indicators
        try
            R.high_freq_high_oscill = double(R.PVCs_per_hour > 30 && R.oscill_rr_median > 0.3);
            R.slow_recovery_high_cv = double(R.halflife_rr30_median > 20 && R.halflife_cv > 0.5);
            R.composite_high_risk = double(R.composite_risk_score > 0.6);
        catch
            R.high_freq_high_oscill = NaN;
            R.slow_recovery_high_cv = NaN;
            R.composite_high_risk = NaN;
        end
        
        % To table
        T = struct2table(R, 'AsArray', true);
    end

    % Read patient vital info (vital: 0=Dead,1=Alive)
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
        % Candidate path
        p1 = fullfile(pwd,'shhs','datasets','shhs-cvd-summary-dataset-0.21.0.csv');
        cand = p1;
        if ~exist(cand,'file')
            d1 = dir(fullfile(pwd,'shhs','datasets','**','shhs-cvd-summary-dataset-*.csv'));
            if ~isempty(d1)
                % Take latest by modification time
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
        % Match
        idx = find(ids == idNum, 1, 'first');
        if isempty(idx), return; end
        % Parse vital
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

    % Summarize covariates (same naming as training: age_s1, gender, race, censdate)
    function demo = getSummaryFieldsForCurrentRecord()
        demo = struct('age_s1', NaN, 'gender', NaN, 'race', NaN, 'censdate', NaN);
        try
            [~, base, ~] = fileparts(state.edfPath);
            tok = regexp(base, '([0-9]+)$', 'tokens', 'once');
            if isempty(tok)
                tok = regexp(base, '([0-9]{6,})', 'tokens', 'once');
            end
            if isempty(tok), return; end
            idStr = tok{1}; idNum = str2double(idStr);
            p1 = fullfile(pwd,'shhs','datasets','shhs-cvd-summary-dataset-0.21.0.csv');
            cand = p1;
            if ~exist(cand,'file')
                d1 = dir(fullfile(pwd,'shhs','datasets','**','shhs-cvd-summary-dataset-*.csv'));
                if ~isempty(d1)
                    [~,mi] = max([d1.datenum]); %#ok<DATNM>
                    cand = fullfile(d1(mi).folder, d1(mi).name);
                else
                    cand = '';
                end
            end
            if isempty(cand) || ~exist(cand,'file'), return; end
            T = readtable(cand);
            vn = lower(T.Properties.VariableNames);
            iId = find(strcmp(vn,'nsrrid'),1); if isempty(iId), iId = find(contains(vn,'nsrrid'),1); end
            if isempty(iId), return; end
            idsCol = T{:,iId};
            if ~isnumeric(idsCol), ids = str2double(string(idsCol)); else, ids = double(idsCol); end
            idx = find(ids == idNum, 1, 'first');
            if isempty(idx), return; end
            % age_s1
            iAge = find(strcmp(vn,'age_s1'),1); if isempty(iAge), iAge = find(contains(vn,'age'),1); end
            if ~isempty(iAge)
                v = T{idx,iAge}; if ~isnumeric(v), v = str2double(string(v)); end
                if isfinite(v), demo.age_s1 = double(v); end
            end
            % gender
            iG = find(strcmp(vn,'gender'),1); if ~isempty(iG)
                v = T{idx,iG}; if ~isnumeric(v), v = str2double(string(v)); end
                if isfinite(v), demo.gender = double(v); end
            end
            % race
            iR = find(strcmp(vn,'race'),1); if ~isempty(iR)
                v = T{idx,iR}; if ~isnumeric(v), v = str2double(string(v)); end
                if isfinite(v), demo.race = double(v); end
            end
            % censdate
            iC = find(strcmp(vn,'censdate'),1); if ~isempty(iC)
                v = T{idx,iC}; if ~isnumeric(v), v = str2double(string(v)); end
                if isfinite(v), demo.censdate = double(v); end
            end
        catch
        end
    end

    % Align feature table with model predictorNames: rename/derive/fill as necessary
    function Trec = alignFeaturesForSCA(Trec, predictorNames)
        % Feature rename mapping (aligned with generate_sca_classifier_v3.m lines 856-886)
        % Map GUI-generated original names to model-expected ECG statistic names
        featureNameMap = containers.Map(...
            {'oscill_pvc_log_ratio', 'oscill_rr_median', 'oscill_x_cv', 'oscill_log', ...
             'recovery_capacity', 'oscill_log_x_pvc', 'tau_rr_median', 'oscill_x_pvc_freq', ...
             'early_late_ratio', 'halflife_cv', 'HRT_TO_abnormal_frac', 'pvc_burden_index', ...
             'PVCs_per_hour', 'pvc_freq_sqrt', 'composite_risk_score'}, ...
            {'RRI_Oscill_Freq_LogRatio', 'RRI_PostPVC_OscillAmp_Median', 'RRI_Oscill_x_RecoveryCV', 'RRI_PostPVC_OscillAmp_Log', ...
             'RRI_Recovery_CapacityIndex', 'RRI_OscillLog_x_PVC_Freq', 'RRI_Recovery_TimeConst_Median', 'RRI_Oscill_x_PVC_Freq', ...
             'RRI_Recovery_EarlyLate_Ratio', 'RRI_Recovery_TimeVariability_CV', 'HRT_TurbOnset_AbnormalFrac', 'PVC_Burden_NormIndex', ...
             'PVC_FreqPerHour', 'PVC_FreqPerHour_Sqrt', 'Clinical_CompositeRisk_Score'});
        
        % Rename GUI-generated features
        oldNames = featureNameMap.keys;
        for i = 1:numel(oldNames)
            oldName = oldNames{i};
            newName = featureNameMap(oldName);
            % If GUI table has old name and model requires new name
            if ismember(oldName, Trec.Properties.VariableNames) && ismember(newName, predictorNames)
                Trec.(newName) = Trec.(oldName);
            end
        end
        
        % Fill missing required features
        for i = 1:numel(predictorNames)
            fn = predictorNames{i};
            if ~ismember(fn, Trec.Properties.VariableNames)
                Trec.(fn) = NaN;
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
    
    % ====== Recovery dynamics helper functions (aligned with generate_sca_classifier_v3.m) ======
    function t50 = local_halflife(dev, t, thr, consecBeats)
        % Compute time when deviation recovers below thr and stays for consecBeats beats
        t50 = NaN;
        if isempty(dev)||isempty(t), return; end
        dev = dev(:); t = t(:);
        hit = dev <= thr;
        run = 0;
        for i = 1:numel(hit)
            if hit(i), run=run+1; else, run=0; end
            if run>=consecBeats
                t50 = t(i - consecBeats + 1); return;
            end
        end
    end
    
    function tau = local_time_constant(dev, t, K)
        % Recovery time constant (exponential fit)
        tau = NaN;
        if isempty(dev)||isempty(t), return; end
        dev = dev(:); t = t(:);
        K = min([K, numel(dev), numel(t)]);
        y = log(max(dev(1:K), eps));
        x = t(1:K);
        if numel(unique(x))<2, return; end
        p = polyfit(x,y,1);
        a = p(1);
        if a < 0, tau = -1/a; end
    end
    
    function A = local_auc(dev, t, normTime)
        % AUC of deviation curve
        A=NaN;
        if numel(dev)<2 || numel(t)<2, A=0; return; end
        dev = dev(:); t = t(:);
        dt = diff(t)./max(normTime,1);
        for i=2:numel(dev)
            if isnan(dev(i)), dev(i)=dev(i-1); end
        end
        if isnan(dev(1)), dev(1)=0; end
        A = sum(0.5*(dev(1:end-1)+dev(2:end)).*dt);
        if ~isfinite(A), A = NaN; end
    end
    
    function oi = local_oscillation_index(dev)
        % Oscillation index: normalized count of sign changes
        oi = NaN;
        if numel(dev)<3, oi=0; return; end
        d1 = diff(dev(:));
        s  = sign(d1); s(s==0)=1;
        oi = sum(abs(diff(s))>0) / max(numel(d1)-1,1);
    end
    
    function cnt = local_recovery_stalls(dev, thr, consecBeats)
        % Number of recovery stalls: re-leaving recovery region after entering
        cnt=0;
        if isempty(dev), return; end
        under = dev(:) <= thr;
        if ~any(under), return; end
        enterIdx = find([false; (~under(1:end-1) & under(2:end))]);
        i=1; validSegments=0;
        while i<=numel(enterIdx)
            startPos = enterIdx(i)+1;
            len=0; j=startPos;
            while j<=numel(under) && under(j)
                len=len+1; j=j+1;
            end
            if len>=consecBeats, validSegments = validSegments+1; end
            i=i+1;
        end
        cnt = max(validSegments-1,0);
    end
    
    function ts_val = local_approx_ts(idxPVC, pvcSample, obsEnd, rGlobalAll, rr_between, isPVCBeat, sqi_vec, numBeats)
        % HRT TS: maximum 5-consecutive slope
        j_after = (idxPVC+1):min(numBeats, idxPVC+20);
        j_after = j_after(:);
        mask_rr_ok = (rGlobalAll(j_after) <= obsEnd) & (rGlobalAll(j_after-1) >= pvcSample) & ...
            ~isPVCBeat(j_after) & ~isPVCBeat(j_after-1) & sqi_vec(j_after) & sqi_vec(j_after-1) & ...
            isfinite(rr_between(j_after));
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
    end
    
    function slope = local_linear_slope(t,x)
        % Linear slope
        slope=NaN;
        t=t(:); x=x(:);
        ok=isfinite(t)&isfinite(x); t=t(ok); x=x(ok);
        if numel(t)<2, return; end
        p=polyfit(t,x,1); slope=p(1);
    end
    
    function s = local_sampen(x,m,r)
        % Sample entropy
        s=NaN;
        x=x(:); x=x(isfinite(x));
        n=numel(x);
        if n < m+2 || ~isfinite(r) || r<=0, return; end
        count_m=0; count_m1=0;
        for i=1:(n-m)
            Xi = x(i:(i+m-1));
            for j=(i+1):(n-m+1)
                Xj = x(j:(j+m-1));
                if max(abs(Xi-Xj))<=r
                    count_m = count_m + 1;
                    if j <= n - m
                        if abs(x(i+m)-x(j+m))<=r
                            count_m1 = count_m1 + 1;
                        end
                    end
                end
            end
        end
        if count_m==0 || count_m1==0, return; end
        s = -log(count_m1 / count_m);
    end
    
    function c = local_lzc_binary(b)
        % Lempel-Ziv complexity (binary)
        c=NaN;
        b=b(:); b=b(isfinite(b));
        if isempty(b), return; end
        b=double(b~=0);
        n=numel(b);
        if n<2, c=0; return; end
        s=char(b+'0').';
        i=1; c_raw=1; l=1; k=1; k_max=1; nlen=length(s);
        while true
            if i+k>nlen || l+k>nlen
                c_raw=c_raw+1; break;
            end
            if s(i+k)==s(l+k)
                k=k+1; if k>k_max, k_max=k; end
            else
                if k>k_max, k_max=k; end
                if k_max==1
                    c_raw=c_raw+1; l=l+1; i=1; k=1; k_max=1;
                else
                    i=i+1;
                    if i==l
                        c_raw=c_raw+1; l=l+k_max; i=1; k=1; k_max=1;
                    else
                        k=1;
                    end
                end
            end
            if l>nlen, break; end
        end
        c = c_raw * (log2(n))/n;
    end
    
    function [SD1, SD2] = local_poincare_sd(rr)
        % Poincaré SD1/SD2
        SD1=NaN; SD2=NaN;
        rr=rr(:); rr=rr(isfinite(rr));
        if numel(rr)<3, return; end
        d=diff(rr);
        sd_rr=std(rr); sd_d=std(d);
        SD1 = sqrt(0.5)*sd_d;
        tmp = 2*sd_rr^2 - 0.5*sd_d^2;
        if tmp>0, SD2 = sqrt(tmp); end
    end
    
    function v = local_safe_max(x)
        x=x(:); x=x(isfinite(x));
        if isempty(x), v=NaN; else, v=max(x); end
    end
    
    function v = local_safe_min(x)
        x=x(:); x=x(isfinite(x));
        if isempty(x), v=NaN; else, v=min(x); end
    end
end
