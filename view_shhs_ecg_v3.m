function view_shhs_ecg_v3()
% File: view_shhs_ecg_v3.m
% Type: Function (GUI tool)
% Description:
%   Integrated browsing and analysis: predict PVCs with filtering/locating, aggregate post-PVC 5-second features per record,
%   and estimate record-level SCA risk using a trained model.
% Usage:
%   Run in MATLAB: view_shhs_ecg_v3
% Dependencies:
%   edfread, ecgFilter, detectAndClassifyHeartbeats, extractHeartbeatFeatures
%   and trained models in the results directory
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26

    % Shared state
    state = struct();
    state.fs = 125;                % default SHHS1 sampling rate
    state.windowSec = 30;          % browsing window length (seconds)
    state.viewStartSec = 0;        % current window start (seconds)
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
    state.clsModel = [];      % beat classifier model
    state.scaModel = [];      % SCA risk model
    state.scaResultText = '';
    state.patientVital = NaN; % 0=Dead,1=Alive

    % Build UI
    ui = buildUI();
    updateUIState('init');

    % Internal: build UI
    function ui = buildUI()
        screenSize = get(0, 'ScreenSize');
        figW = min(1200, screenSize(3) * 0.9);
        figH = min(720, screenSize(4) * 0.85);
        ui.fig = figure('Name','SHHS ECG PVC & SCA Risk Analysis', 'NumberTitle','off', ...
            'MenuBar','none', 'ToolBar','none', 'Color','w', 'Units','pixels', ...
            'Position',[100 100 figW figH]);

        % Layout: top main axis, bottom control area
        ctrlH = 230;
        ui.ax = axes('Parent', ui.fig, 'Units','pixels', 'Position', [60 ctrlH+30 figW-100 figH-ctrlH-80]);
        set(ui.ax,'Color','w');
        grid(ui.ax,'on'); hold(ui.ax,'on');
        xlabel(ui.ax,'Time (s)'); ylabel(ui.ax,'ECG (a.u.)');

        % Control area (using normalized layout)
        panel = uipanel('Parent', ui.fig, 'Units','pixels', 'Position',[10 10 figW-20 ctrlH-20], 'BorderType','none', 'BackgroundColor','w');

        % Row 1: import and basic info
        y1 = ctrlH - 50;
        ui.btnLoad = uicontrol(panel,'Style','pushbutton','String','Select EDF…','Units','pixels', ...
            'Position',[10 y1 120 28],'Callback',@onLoadEDF,'BackgroundColor','w');
        ui.textFile = uicontrol(panel,'Style','text','String','File: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y1 460 22],'BackgroundColor','w');
        ui.textSize = uicontrol(panel,'Style','text','String','Size: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[610 y1 160 22],'BackgroundColor','w');
        ui.textFs   = uicontrol(panel,'Style','text','String','fs: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[780 y1 120 22],'BackgroundColor','w');
        ui.textDur  = uicontrol(panel,'Style','text','String','Duration: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[910 y1 200 22],'BackgroundColor','w');

        % Row 2: browsing and positioning
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

        % Row 4: navigation and filtering
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
        ui.textSCA = uicontrol(panel,'Style','text','String','SCA risk: -','HorizontalAlignment','left', ...
            'Units','pixels','Position',[140 y5 960 22],'BackgroundColor','w');

        % Pre-created plot layer handles
        ui.hECG = plot(ui.ax, NaN, NaN, 'b-','LineWidth',1); % ECG main line
        ui.hR_PVC = plot(ui.ax, NaN, NaN, 'rv','MarkerFaceColor','r','MarkerSize',6,'LineStyle','none'); % PVC R
        ui.hR_Other = plot(ui.ax, NaN, NaN, 'k^','MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',5,'LineStyle','none'); % Other R
        ui.hSel = plot(ui.ax, NaN, NaN, 'go','MarkerSize',10,'LineWidth',2,'LineStyle','none'); % selected highlight
        ui.textGroup = gobjects(0); % text labels

        % Enable built-in zoom/pan and refresh accordingly
        ui.zoomObj = zoom(ui.fig); ui.zoomObj.Motion = 'both'; ui.zoomObj.Enable = 'on';
        ui.panObj = pan(ui.fig); ui.panObj.Enable = 'on';
        ui.zoomObj.ActionPostCallback = @onZoomPan;
        ui.panObj.ActionPostCallback = @onZoomPan;
        set(ui.fig,'WindowScrollWheelFcn', @onScrollZoom);
    end

    % Mouse wheel zoom (around cursor position)
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
        % Button availability
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
            set(ui.textFs,'String','fs: -');
            set(ui.textDur,'String','Duration: -');
        else
            [~,nm,ext] = fileparts(state.edfPath);
            set(ui.textFile,'String',sprintf('File: %s%s', nm, ext));
            if isfinite(state.fileSizeBytes)
                set(ui.textSize,'String',sprintf('Size: %.2f MB', state.fileSizeBytes/1024/1024));
            else
                set(ui.textSize,'String','Size: -');
            end
            set(ui.textFs,'String',sprintf('fs: %d Hz', state.fs));
            durSec = state.N / max(state.fs,1);
            set(ui.textDur,'String',sprintf('Duration: %.1f s', durSec));
        end

        % Stats
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
            set(ui.textSCA,'String','SCA risk: -');
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
                errordlg(['ECG channel type not supported: ' class(ecgCol)], 'Channel Error');
                return;
            end
            state.ecg = double(ecg);
            state.N = numel(state.ecg);
            state.time = (0:state.N-1)'/state.fs;

            % Dual filtering: display (2), detection (3); SHHS1 already notched at 60 Hz
            try
                [state.ecgVisFiltered, ~] = ecgFilter(state.ecg, state.fs, 2, 0);
            catch ME
                warndlg(['Display filtering failed (using raw signal): ' ME.message], 'Filtering Warning');
                state.ecgVisFiltered = state.ecg;
            end
            try
                [state.ecgDetFiltered, ~] = ecgFilter(state.ecg, state.fs, 2, 0);
            catch ME
                warndlg(['Detection filtering failed (fallback to display signal): ' ME.message], 'Filtering Warning');
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

    % Set window start and refresh
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
            errordlg('Please import an EDF first.','Notice');
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
                errordlg('No valid beats obtained.','Detection Failed');
                return;
            end
            state.beatInfo = beatInfo;
            state.stats = stats;

            % Feature extraction
            [featureTable, ~] = extractHeartbeatFeatures(beatInfo, state.fs);
            if isempty(featureTable) || height(featureTable) == 0
                errordlg('Feature extraction is empty.','Error');
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

            % Select required variables and impute missing
            if isfield(trainedClassifier,'RequiredVariables') && ~isempty(trainedClassifier.RequiredVariables)
                reqVars = trainedClassifier.RequiredVariables;
            else
                reqVars = setdiff(featureTable.Properties.VariableNames, {'BeatType'});
            end
            miss = setdiff(reqVars, featureTable.Properties.VariableNames);
            if ~isempty(miss)
                errordlg(['Features missing required variables: ' strjoin(miss, ', ')],'Missing Variables');
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

            % Prediction
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

            % Filter and select
            rebuildFilteredIndices();
            if ~isempty(state.filteredBeatIndices)
                state.selectedBeatGlobalIdx = state.filteredBeatIndices(1);
                centerOnSelectedBeat();
            end

            updateUIState('predicted_pvc');
        catch ME
            errordlg(['PVC detection pipeline failed: ' ME.message], 'Error');
        end
    end

    % Predict SCA risk (record-level aggregated features)
    function onPredictSCA(~,~)
        if isempty(state.rGlobalAll) || ~any(state.pvcMask)
            warndlg('No PVC beats available to predict SCA risk.','Notice');
            return;
        end
        try
            % Load SCA model
            modelFile = fullfile(pwd,'results','sca_classifier_post_ectopic.mat');
            if ~exist(modelFile,'file')
                errordlg(['SCA model not found: ' modelFile],'Model Missing');
                return;
            end
            M = load(modelFile);
            scaModel = M.trainedClassifier;
            predictorNames = M.predictorNames;
            state.scaModel = scaModel;

            % Compute record-level aggregated features (consistent with training)
            pvcRidx = state.rGlobalAll(state.pvcMask);
            Trec = computePostPVCFeatureTable(pvcRidx, 5.0); % Now single row aggregated table
            if isempty(Trec) || height(Trec)~=1
                errordlg('Record-level aggregated features not obtained.','Error');
                return;
            end

            % Select predictors and impute
            miss = setdiff(predictorNames, Trec.Properties.VariableNames);
            if ~isempty(miss)
                errordlg(['Missing SCA model required variables: ' strjoin(miss, ', ')],'Missing Variables');
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

            % Predict death risk score (scores = P(y==0))
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

            % Query patient survival info
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
            % state.scaResultText = sprintf('SCA risk: %s | PVCs=%d, PVC density=%.2f /h | Death risk=%.3f, Thr=%.2f | Actual=%s', ...
            %     ternary(hasRisk,'At risk','No risk'), pvcCount, pvcPerHour, scoreDead, thr, vitalStr);
            state.scaResultText = sprintf('SCA risk: %s | Actual=%s', ...
            ternary(hasRisk,'At risk','No risk'), vitalStr);

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
            % Use detection-filtered indices for position, amplitude from display-filtered signal
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

    % Interpolate selected point's y (bounds-safe)
    function yv = interpY(iStart, yseg, globalIdx)
        % globalIdx is index in full signal; convert to segment index
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

    % Current selection position within filtered sequence
    function pos = findCurrentFilteredPos()
        pos = NaN;
        if isempty(state.filteredBeatIndices), return; end
        pos = find(state.filteredBeatIndices == state.selectedBeatGlobalIdx, 1, 'first');
        if isempty(pos)
            % If not in list, choose nearest
            dif = abs(state.filteredBeatIndices - state.selectedBeatGlobalIdx);
            [~,mi] = min(dif);
            pos = mi;
        end
    end

    % Step selection left/right
    function stepSelection(step)
        if isempty(state.filteredBeatIndices), return; end
        curPos = findCurrentFilteredPos();
        if ~isfinite(curPos), curPos = 1; end
        newPos = max(1, min(numel(state.filteredBeatIndices), curPos + step));
        state.selectedBeatGlobalIdx = state.filteredBeatIndices(newPos);
        centerOnSelectedBeat();
        updateUIState('step_sel');
    end

    % Filter change
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

    % Compute record-level 5-second post-PVC feature aggregation (single row)
    function T = computePostPVCFeatureTable(pvcIndices, postWindowSec)
        if isempty(pvcIndices)
            T = table(); return;
        end
        fs = state.fs; N = state.N;
        rGlobalAll = state.rGlobalAll(:);
        % Map to beat indices (nearest neighbor)
        pvcIndicesSorted = sort(pvcIndices(:));
        idxPVC_all_sorted = round(interp1(rGlobalAll, 1:numel(rGlobalAll), pvcIndicesSorted, 'nearest', 'extrap'));
        idxPVC_all_sorted = max(1, min(numel(rGlobalAll), idxPVC_all_sorted));

        % RR vectors
        numBeatsRec = numel(rGlobalAll);
        rr_pre_vec = nan(numBeatsRec,1); rr_post1_vec = nan(numBeatsRec,1); rr_post2_vec = nan(numBeatsRec,1);
        if numBeatsRec >= 2
            rr_pre_vec(2:end) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
            rr_post1_vec(1:end-1) = (rGlobalAll(2:end) - rGlobalAll(1:end-1)) / fs;
        end
        if numBeatsRec >= 3
            rr_post2_vec(1:end-2) = (rGlobalAll(3:end) - rGlobalAll(1:end-2)) / fs;
        end

        % QRS duration and R amplitude
        qrs_dur_vec = nan(numBeatsRec,1);
        r_amp_vec = nan(numBeatsRec,1);
        for ii = 1:numBeatsRec
            b = state.beatInfo(ii);
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

        % Window bounds
        winLen = round(postWindowSec * fs);
        winStartVec = pvcIndicesSorted;
        winEndVec = min(N, winStartVec + winLen);
        winDurSecVec = (winEndVec - winStartVec) / fs;

        % Count beats within 5 s window (two pointers)
        counts = zeros(numel(pvcIndicesSorted),1);
        left = 1; right = 0; numPVC = numel(pvcIndicesSorted);
        for kk = 1:numPVC
            a = winStartVec(kk); b = winEndVec(kk);
            while right < numBeatsRec && rGlobalAll(right+1) <= b
                right = right + 1;
            end
            while left <= numBeatsRec && rGlobalAll(left) < a
                left = left + 1;
            end
            if right >= left
                counts(kk) = right - left + 1;
            else
                counts(kk) = 0;
            end
        end

        % Window-level features
        RR_Pre = rr_pre_vec(idxPVC_all_sorted);
        RR_Post1 = rr_post1_vec(idxPVC_all_sorted);
        RR_Post2 = rr_post2_vec(idxPVC_all_sorted);
        HRT_TO = (RR_Post1 - RR_Pre) ./ (RR_Pre + eps);
        HRT_TS = (RR_Post2 - RR_Pre) / 2;
        QRS_Dur_PVC = qrs_dur_vec(idxPVC_all_sorted);
        R_Amp_PVC = r_amp_vec(idxPVC_all_sorted);
        Beats_in_5s = counts;
        HR_Pre = 60 ./ RR_Pre;
        HR_Post1 = 60 ./ RR_Post1;
        HR_5s = (counts ./ max(winDurSecVec, eps)) * 60;

        % Derived vectors
        HR_Accel = HR_Post1 - HR_Pre;
        CompRatio = RR_Post1 ./ (RR_Pre + eps);
        pvcIntervals = diff(pvcIndicesSorted) / fs;

        % Record-level stats
        recordDurationHr = (N / fs) / 3600;
        pvcPerHour = numPVC / max(recordDurationHr, eps);

        % Assemble aggregation struct
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
        S = agg_add_stats(S, HR_Accel,  'HR_Accel');
        S = agg_add_stats(S, CompRatio, 'CompRatio');
        S = agg_add_stats(S, pvcIntervals, 'PVC_Interval');

        % Proportions and HRV-derived metrics
        mask_to = isfinite(HRT_TO);
        if any(mask_to)
            S.HRT_TO_neg_frac = mean(double(HRT_TO(mask_to) < 0));
        else
            S.HRT_TO_neg_frac = NaN;
        end
        mask_qrs = isfinite(QRS_Dur_PVC);
        if any(mask_qrs)
            S.QRS_Prolonged_frac = mean(double(QRS_Dur_PVC(mask_qrs) > 0.12));
        else
            S.QRS_Prolonged_frac = NaN;
        end
        S.RR_Pre_CV   = agg_cv(RR_Pre);
        S.RR_Post1_CV = agg_cv(RR_Post1);
        S.RR_Pre_RMSSD   = agg_rmssd(RR_Pre);
        S.RR_Post1_RMSSD = agg_rmssd(RR_Post1);

        T = struct2table(S);
    end

    % Read patient's actual survival information (vital: 0=Dead,1=Alive)
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
                % Take the latest modification time
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


