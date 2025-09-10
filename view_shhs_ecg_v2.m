function view_shhs_ecg_v2()
% =================================================================================================
% File: view_shhs_ecg_v2.m
% Overview:
%   SHHS ECG (QRS/T) visualization viewer, edition v2, with interactive controls.
%   Dual synchronized axes: top shows filtered ECG for detection with Q/R/S/T markers; bottom shows raw ECG and (if available) CSV R points with types.
%
% Responsibilities:
%   1. Load a single SHHS EDF file (user interactive selection) and extract the ECG channel (match ECG/EKG, supports cell array or numeric column).
%   2. Preprocess with ecgFilter (mode 2: IIR/enhanced pipeline) to produce the detection signal (aligned with the top display).
%   3. Call detectAndClassifyHeartbeats to locate beats, convert relative indices to absolute sample indices; supports Q/R/S/T.
%   4. Attempt to parse corresponding rpoints annotations (multiple filename candidates -rpoint(s).csv, recursive search) and map by seconds or by sampling rate to current fs.
%   5. Interaction: synchronized zoom/pan, time jump (seconds or hh:mm:ss), adaptive decimation to improve plotting performance.
%   6. Info bar shows sampling rate, annotated PVC/Other counts, and detected R peak count.
%
% I/O:
%   External call: no input parameters; run directly to launch the GUI.
%   User input: menu "Open EDF..."; optional time jump edit; mouse wheel/drag for zoom/pan.
%   Data source: SHHS EDF plus optional CSV in annotations-rpoints folder (fields seconds / Type / rpoint).
%   Return: none (interactive figure output only).
%
% Key implementation notes:
%   - Centralized state struct S holds data/handles, passed across callbacks via guidata.
%   - Adaptive sampling: decimate by dynamic step size based on window length and displayMaxPoints to balance performance and clarity.
%   - Marker throttling: limit the number of Q/R/S/T and bottom R/type markers via displayMaxMarkers to prevent GUI slowdowns.
%   - Time jump: parse seconds or hh:mm:ss -> center around the time with a half-window; auto-extend to ensure minimal visible width.
%   - Annotation parsing: prefer seconds; if missing, use rpoint indices with fs mapping; fall back to default fs if sampling rate is missing.
%   - Q/R/S/T index sanitization: sanitizeIndices removes out-of-range/NaN and ensures integers.
%
% Robustness and safeguards:
%   - EDF read failure shows a dialog and aborts without breaking existing state.
%   - Channel search uses multiple strategies (exact ECG match and fuzzy contains ecg/ekg).
%   - Accepts ECG as cell array or numeric vector; attempt str2num conversion for textual numbers.
%   - Beat detection wrapped with try/catch; log warnings only and keep the GUI responsive.
%   - CSV multiple path candidates plus recursive search; lowercase variable-name matching avoids case issues.
%   - Edge cases: empty data, empty window range, or constant signals trigger y-axis expansion for visibility.
%   - Limit text label count (maxLabels) to avoid rendering degradation.
%
% Performance notes:
%   - Plot decimation (limit points and markers).
%   - Use drawnow limitrate to reduce refresh frequency during rapid zooming.
%   - Filter and redraw only markers/labels within the visible window.
%
% Suitable/Not suitable:
%   Suitable: manual QC of a single file, debugging beat detection/annotation quality, quick validation of QRS/T localization.
%   Not suitable: batch statistics or offline processing (should use scripted pipelines).
%
% Changelog:
%   2025-08-30 Initial standardized header comment (Chinese version); functionality unchanged.
%   2025-09-09 Full English translation of comments and UI strings; functionality unchanged.
%
% Usage:
%   >> view_shhs_ecg_v2
%   After launching the GUI, open an EDF via the menu. Use the jump box to enter "3600" or "01:00:00" for quick navigation.
%
% Dependencies:
%   ecgFilter.m, detectAndClassifyHeartbeats.m, and their feature/classification support files.
% =================================================================================================
	% State
	S = struct();
	S.fs_default = 125;
	S.fs = S.fs_default;
	S.ecgDispFull = [];   % Display/detection waveform (filtered, shown on top)
	S.ecgRaw = [];        % Raw, unfiltered ECG (shown on bottom)
	S.tFull = [];
	S.qIdxFull = []; S.rIdxFull = []; S.sIdxFull = []; S.tIdxFull = [];
	S.qTimesFull = []; S.rTimesFull = []; S.sTimesFull = []; S.tTimesFull = [];
	S.rSecsFull = []; S.rTypesFull = {};
	S.edfPath = ''; S.annPath = '';
	S.displayMaxPoints = 20000;
	S.displayMaxMarkers = 5000;
	S.maxLabels = 300;
	S.labelOffsetFrac = [0.02 0.01];
	S.jumpWindowSec = 10;

	% UI
	S.fig = figure('Name','SHHS ECG QRS Viewer', 'NumberTitle','off', 'Color','w', ...
		'Units','normalized', 'Position',[0.08 0.08 0.84 0.84]);
	mFile = uimenu(S.fig, 'Label','File');
	uimenu(mFile, 'Label','Open EDF...', 'Accelerator','O', 'Callback', @(src,evt)onOpenEdf(src));
	uimenu(mFile, 'Label','Clear', 'Callback', @(src,evt)onClear(src));
	uimenu(mFile, 'Label','Exit', 'Callback', @(~,~)close(S.fig));

	% Navigation menu: go to time
	mNav = uimenu(S.fig, 'Label','Navigate');
	uimenu(mNav, 'Label','Go to time...', 'Accelerator','G', 'Callback', @(src,evt)onGotoPrompt(src));

	% Create a panel to host two axes, reserving bottom space for the info bar
	S.axPanel = uipanel('Parent', S.fig, 'Units','normalized', 'Position',[0.05 0.10 0.90 0.85], ...
		'BorderType','none', 'BackgroundColor','w');
	S.tiled = tiledlayout(S.axPanel, 2, 1, 'TileSpacing','compact', 'Padding','compact');
	S.axTop = nexttile(S.tiled, 1); hold(S.axTop,'on'); grid(S.axTop,'on'); box(S.axTop,'on');
	title(S.axTop,'ECG with QRS Detection Viewer'); xlabel(S.axTop,'Time (s)'); ylabel(S.axTop,'Amplitude');
	S.axBottom = nexttile(S.tiled, 2); hold(S.axBottom,'on'); grid(S.axBottom,'on'); box(S.axBottom,'on');
	title(S.axBottom,'Raw ECG (unfiltered) with R annotations'); xlabel(S.axBottom,'Time (s)'); ylabel(S.axBottom,'Amplitude');

	% Lines and markers
	S.hTop = plot(S.axTop, NaN, NaN, 'Color',[0.2 0.2 0.2], 'DisplayName','ECG');
	S.hBottom = plot(S.axBottom, NaN, NaN, 'b-');
	S.hQ = plot(S.axTop, NaN, NaN, 'g<', 'MarkerSize',5, 'LineWidth',1.0, 'LineStyle','none', 'DisplayName','Q');
	S.hR = plot(S.axTop, NaN, NaN, 'r^', 'MarkerSize',6, 'LineWidth',1.2, 'LineStyle','none', 'DisplayName','R');
	S.hS = plot(S.axTop, NaN, NaN, 'b>', 'MarkerSize',5, 'LineWidth',1.0, 'LineStyle','none', 'DisplayName','S');
	S.hT = plot(S.axTop, NaN, NaN, 'mv', 'MarkerSize',5, 'LineWidth',1.0, 'LineStyle','none', 'DisplayName','T');
	S.hRBottom = plot(S.axBottom, NaN, NaN, 'rv', 'MarkerSize',6, 'LineWidth',1.0);
	S.hTxtBottom = gobjects(0);

	% Axes linkage and interactions
	linkaxes([S.axTop, S.axBottom], 'x');
	S.zoomObj = zoom(S.fig); S.zoomObj.Motion='both'; S.zoomObj.Enable='on';
	S.panObj = pan(S.fig); S.panObj.Enable='on';
	S.zoomObj.ActionPostCallback = @(obj,evd)onZoomPan(evd);
	S.panObj.ActionPostCallback = @(obj,evd)onZoomPan(evd);

	guidata(S.fig, S);

	% Bottom info bar (file name, sampling rate, annotation and detection stats)
	S = guidata(S.fig);
	S.hInfoText = uicontrol(S.fig, 'Style','text', 'Units','normalized', ...
		'Position',[0.02 0.01 0.62 0.06], 'BackgroundColor','w', ...
		'HorizontalAlignment','left', 'FontSize',10, 'String','');
	% Bottom quick-jump controls
	S.hJumpLabel = uicontrol(S.fig, 'Style','text', 'Units','normalized', ...
		'Position',[0.64 0.01 0.10 0.06], 'BackgroundColor','w', ...
		'HorizontalAlignment','right', 'String','Go to time:');
	S.hJumpEdit = uicontrol(S.fig, 'Style','edit', 'Units','normalized', ...
		'Position',[0.75 0.02 0.10 0.04], 'String','00:00:10', ...
		'TooltipString','Enter seconds or hh:mm:ss, press Enter or click "Go"', ...
		'KeyPressFcn', @(obj,evt)onJumpEditKey(obj,evt));
	S.hJumpBtn = uicontrol(S.fig, 'Style','pushbutton', 'Units','normalized', ...
		'Position',[0.86 0.02 0.10 0.04], 'String','Go', ...
		'Callback', @(src,evt)onJumpClick(src));
	guidata(S.fig, S);

	function onOpenEdf(src)
		fig = ancestor(src,'figure'); S = guidata(fig);
		startDir = fullfile(pwd, 'shhs','polysomnography','edfs');
		if ~isfolder(startDir), startDir = pwd; end
		[fn, fd] = uigetfile({'*.edf','EDF Files (*.edf)'}, 'Select EDF file', startDir);
		if isequal(fn,0), return; end
		edfPath = fullfile(fd, fn);
		try
			TT = edfread(edfPath);
		catch ME
			errordlg(sprintf('edfread failed:\n%s', ME.message), 'Read Error');
			return;
		end
		% Locate ECG channel
		vn = TT.Properties.VariableNames;
		iEcg = find(strcmp(vn,'ECG'),1); if isempty(iEcg), vl = lower(vn); iEcg = find(contains(vl,'ecg')|contains(vl,'ekg'),1,'first'); end
		if isempty(iEcg), errordlg('ECG channel not found.','Missing Channel'); return; end
		ecgVar = vn{iEcg};
		ecgCol = TT.(ecgVar);
		if iscell(ecgCol)
			try
				ecg = vertcat(ecgCol{:});
			catch
				ecg = [];
				for kk = 1:numel(ecgCol)
					v = ecgCol{kk}; if isstring(v)||ischar(v), v = str2num(v); end %#ok<ST2NM>
					ecg = [ecg; v(:)]; %#ok<AGROW>
				end
			end
		elseif isnumeric(ecgCol)
			ecg = ecgCol(:);
		else
			errordlg(sprintf('Unsupported ECG channel type: %s', class(ecgCol)), 'Unsupported Type'); return;
		end
		ecg = double(ecg);
		fs = S.fs_default; t = (0:numel(ecg)-1)'/fs;
		% Display/detection signal: top uses filtered signal; bottom shows raw signal
		try
			ecgFiltered = ecgFilter(ecg, fs, 2, 0);
		catch ME
			errordlg(sprintf('Filtering failed:\n%s', ME.message), 'Filtering Error'); return;
		end
		detSig = ecgFiltered; % Top display and detection signal are the same

		% Run QRS detection (reference main_v3 flow)
		ATRTIMED = []; ANNOTD = {};
		qIdxAbs = []; %#ok<NASGU>
		rIdxAbs = []; %#ok<NASGU>
		sIdxAbs = []; %#ok<NASGU>
		try
			[~, beatInfo, ~] = detectAndClassifyHeartbeats(detSig, ATRTIMED, ANNOTD, fs);
			if ~isempty(beatInfo)
				% Compatible with struct array fields
				nb = numel(beatInfo);
				qIdxAbs = nan(nb,1); rIdxAbs = nan(nb,1); sIdxAbs = nan(nb,1); tIdxAbs = nan(nb,1);
				for ii = 1:nb
					seg0 = safeField(beatInfo(ii),'segmentStartIndex',0);
					qRel = safeField(beatInfo(ii),'qIndex',NaN);
					rRel = safeField(beatInfo(ii),'rIndex',NaN);
					sRel = safeField(beatInfo(ii),'sIndex',NaN);
					tRel = safeField(beatInfo(ii),'tIndex',NaN);
					if ~isnan(qRel) && qRel>0, qIdxAbs(ii) = seg0 + qRel - 1; end
					if ~isnan(rRel) && rRel>0, rIdxAbs(ii) = seg0 + rRel - 1; end
					if ~isnan(sRel) && sRel>0, sIdxAbs(ii) = seg0 + sRel - 1; end
					if ~isnan(tRel) && tRel>0, tIdxAbs(ii) = seg0 + tRel - 1; end
				end
				qIdxAbs = sanitizeIndices(qIdxAbs, numel(detSig));
				rIdxAbs = sanitizeIndices(rIdxAbs, numel(detSig));
				sIdxAbs = sanitizeIndices(sIdxAbs, numel(detSig));
				tIdxAbs = sanitizeIndices(tIdxAbs, numel(detSig));
			else
				qIdxAbs = []; rIdxAbs = []; sIdxAbs = []; tIdxAbs = [];
			end
		catch ME
			warning(ME.identifier, '%s', ME.message);
			qIdxAbs = []; rIdxAbs = []; sIdxAbs = []; tIdxAbs = [];
		end

		% CSV annotations (for bottom R and types)
		[rIdxAnn, rTypesAnn, annFs, annPath, rSecs] = tryLoadAnnotation(edfPath);
		if ~isempty(rSecs)
			rSecsFull = rSecs(:)';
		elseif ~isempty(rIdxAnn)
			rIdxAnn = mapIndicesToFs(rIdxAnn, annFs, fs, numel(detSig));
			rSecsFull = (rIdxAnn-1)/fs;
		else
			rSecsFull = [];
			rTypesAnn = {};
			annFs = []; %#ok<NASGU>
			annPath = '';
		end

		% Save state
		S.fs = fs; S.tFull = t; S.ecgDispFull = detSig; S.ecgRaw = ecg; S.edfPath = edfPath; S.annPath = annPath;
		S.qIdxFull = qIdxAbs(:)'; S.rIdxFull = rIdxAbs(:)'; S.sIdxFull = sIdxAbs(:)'; S.tIdxFull = tIdxAbs(:)';
		S.qTimesFull = (S.qIdxFull-1)/fs; S.rTimesFull = (S.rIdxFull-1)/fs; S.sTimesFull = (S.sIdxFull-1)/fs; S.tTimesFull = (S.tIdxFull-1)/fs;
		S.rSecsFull = rSecsFull; S.rTypesFull = rTypesAnn;
		guidata(S.fig, S);

		% Update bottom info bar
		S = guidata(S.fig);
		[~, edfBaseName, edfExt] = fileparts(edfPath);
		numDetectedR = numel(S.rIdxFull);
		numPVC = 0; numOther = 0; numTotal = 0;
		if ~isempty(S.rTypesFull)
			numPVC = sum(strcmpi(S.rTypesFull, 'PVC'));
			numOther = sum(strcmpi(S.rTypesFull, 'Other'));
			numTotal = numel(S.rTypesFull);
		elseif ~isempty(S.rSecsFull)
			numTotal = numel(S.rSecsFull);
		end
		infoStr = sprintf('File: %s | Sampling rate: %g Hz | Annotations: PVC=%d, Other=%d, Total beats=%d | Detected R peaks=%d', ...
			[edfBaseName edfExt], fs, numPVC, numOther, numTotal, numDetectedR);
		if isfield(S,'hInfoText') && ishghandle(S.hInfoText)
			set(S.hInfoText, 'String', infoStr);
		end

		% Initial 30-second window
		xStart = t(1); xEnd = min(t(end), xStart + 30);
		xlim(S.axTop, [xStart xEnd]); xlim(S.axBottom, [xStart xEnd]);
		% Top legend
		try
			legend(S.axTop, 'show', 'Location','northoutside', 'Orientation','horizontal');
		catch
		end
		refreshView(S.fig);
	end

	function onClear(src)
		fig = ancestor(src,'figure'); S = guidata(fig);
		set(S.hTop, 'XData',NaN,'YData',NaN); set(S.hBottom,'XData',NaN,'YData',NaN);
		set(S.hQ,'XData',NaN,'YData',NaN); set(S.hR,'XData',NaN,'YData',NaN); set(S.hS,'XData',NaN,'YData',NaN); if isfield(S,'hT'), set(S.hT,'XData',NaN,'YData',NaN); end
		set(S.hRBottom,'XData',NaN,'YData',NaN); deleteValidHandles(S.hTxtBottom); S.hTxtBottom = gobjects(0);
		S.fs=S.fs_default; S.ecgDispFull=[]; S.tFull=[]; S.qIdxFull=[]; S.rIdxFull=[]; S.sIdxFull=[]; S.tIdxFull=[];
		S.qTimesFull=[]; S.rTimesFull=[]; S.sTimesFull=[]; S.tTimesFull=[]; S.rSecsFull=[]; S.rTypesFull={}; S.edfPath=''; S.annPath='';
		guidata(S.fig, S);
	end

	function onZoomPan(evd)
		fig = ancestor(evd.Axes,'figure'); S = guidata(fig);
		if isempty(S) || ~isfield(S,'axTop') || ~ishandle(S.axTop), return; end
		axSrc = evd.Axes; axDst = S.axBottom; if axSrc==S.axBottom, axDst=S.axTop; end
		refreshView(fig);
		try
			set(axDst,'YLim', get(axSrc,'YLim'));
		catch
		end
	end

	function refreshView(fig)
		S = guidata(fig); if isempty(S) || isempty(S.tFull), return; end
		t = S.tFull; fs = S.fs; n = numel(t);
		xl = get(S.axTop,'XLim'); x1 = max(t(1), xl(1)); x2 = min(t(end), xl(2)); if x2<=x1, x2=min(t(end), x1+1/fs); end
		iStart = max(1, floor(x1*fs)+1); iEnd = min(n, ceil(x2*fs)+1); if iEnd<=iStart, iEnd=min(n, iStart+1); end
		segLen = iEnd-iStart+1; step = max(1, ceil(segLen / S.displayMaxPoints)); idx = iStart:step:iEnd;
		xp = t(idx); ypTop = S.ecgDispFull(idx); ypBottom = S.ecgRaw(idx);
		set(S.hTop,'XData',xp,'YData',ypTop); set(S.hBottom,'XData',xp,'YData',ypBottom);
		% QRS markers on top axis
		if ~isempty(S.rIdxFull)
			mask = S.rIdxFull>=iStart & S.rIdxFull<=iEnd; rIdxV=S.rIdxFull(mask); rTimesV = S.rTimesFull(mask);
			mstep = max(1, ceil(nnz(mask) / S.displayMaxMarkers)); rIdxV = rIdxV(1:mstep:end); rTimesV = rTimesV(1:mstep:end);
			yR = S.ecgDispFull(rIdxV); set(S.hR,'XData',rTimesV,'YData',yR);
		else
			set(S.hR,'XData',NaN,'YData',NaN);
		end
		if ~isempty(S.qIdxFull)
			maskQ = S.qIdxFull>=iStart & S.qIdxFull<=iEnd; qIdxV=S.qIdxFull(maskQ); qTimesV=S.qTimesFull(maskQ);
			mstepQ = max(1, ceil(nnz(maskQ) / S.displayMaxMarkers)); qIdxV=qIdxV(1:mstepQ:end); qTimesV=qTimesV(1:mstepQ:end);
			yQ = S.ecgDispFull(qIdxV); set(S.hQ,'XData',qTimesV,'YData',yQ);
		else
			set(S.hQ,'XData',NaN,'YData',NaN);
		end
		if ~isempty(S.sIdxFull)
			maskS = S.sIdxFull>=iStart & S.sIdxFull<=iEnd; sIdxV=S.sIdxFull(maskS); sTimesV=S.sTimesFull(maskS);
			mstepS = max(1, ceil(nnz(maskS) / S.displayMaxMarkers)); sIdxV=sIdxV(1:mstepS:end); sTimesV=sTimesV(1:mstepS:end);
			yS = S.ecgDispFull(sIdxV); set(S.hS,'XData',sTimesV,'YData',yS);
		else
			set(S.hS,'XData',NaN,'YData',NaN);
		end
		% T-wave markers on top axis
		if ~isempty(S.tIdxFull)
			maskT = S.tIdxFull>=iStart & S.tIdxFull<=iEnd; tIdxV=S.tIdxFull(maskT); tTimesV=S.tTimesFull(maskT);
			mstepT = max(1, ceil(nnz(maskT) / S.displayMaxMarkers)); tIdxV=tIdxV(1:mstepT:end); tTimesV=tTimesV(1:mstepT:end);
			yT = S.ecgDispFull(tIdxV); set(S.hT,'XData',tTimesV,'YData',yT);
		else
			set(S.hT,'XData',NaN,'YData',NaN);
		end
		% Bottom R and types (from CSV, if present)
		if ~isempty(S.rSecsFull)
			maskRB = S.rSecsFull>=x1 & S.rSecsFull<=x2;
			rTimesB = S.rSecsFull(maskRB); rTypesB = S.rTypesFull(maskRB);
			mstepB = max(1, ceil(numel(rTimesB) / S.displayMaxMarkers)); rTimesB = rTimesB(1:mstepB:end); rTypesB = rTypesB(1:mstepB:end);
			yRB = interp1(t, S.ecgRaw, rTimesB, 'nearest','extrap');
			set(S.hRBottom,'XData', rTimesB, 'YData', yRB);
			deleteValidHandles(S.hTxtBottom);
			if numel(rTimesB) <= S.maxLabels
				S.hTxtBottom = drawTypeLabels(S.axBottom, rTimesB, yRB, rTypesB, S.labelOffsetFrac);
			else
				S.hTxtBottom = gobjects(0);
			end
		else
			set(S.hRBottom,'XData',NaN,'YData',NaN); deleteValidHandles(S.hTxtBottom); S.hTxtBottom = gobjects(0);
		end
		drawnow limitrate; guidata(fig, S);
	end

	function onGotoPrompt(src)
		fig = ancestor(src,'figure'); S = guidata(fig);
		if isempty(S) || isempty(S.tFull)
			errordlg('Please open an EDF file first.','No Data'); return;
		end
		try
			defStr = get(S.hJumpEdit,'String');
		catch
			defStr = '';
		end
		answer = inputdlg('Enter time (seconds or hh:mm:ss):', 'Go to time', 1, {defStr});
		if isempty(answer), return; end
		try
			set(S.hJumpEdit,'String',answer{1});
		catch
		end
		guidata(fig, S);
		doJumpTime(fig, answer{1});
	end

	function onJumpEditKey(obj, evt)
		if isstruct(evt) && isfield(evt,'Key') && (strcmpi(evt.Key,'return') || strcmpi(evt.Key,'enter'))
			onJumpClick(obj);
		end
	end

	function onJumpClick(src)
		fig = ancestor(src,'figure'); S = guidata(fig);
		if isempty(S), return; end
		try
			tstr = get(S.hJumpEdit,'String');
		catch
			return;
		end
		if isempty(tstr), return; end
		doJumpTime(fig, tstr);
	end

	function secs = parseTimeStr(str)
		str = strtrim(string(str));
		if strlength(str)==0, secs = NaN; return; end
		val = str2double(str);
		if ~isnan(val)
			secs = double(val);
			return;
		end
		partsStr = split(str, ':');
		if isscalar(partsStr)
			secs = NaN; return;
		end
		% Allow fractional seconds in the last segment
		lastSec = str2double(partsStr(end));
		if isnan(lastSec), secs = NaN; return; end
		if numel(partsStr)==2
			mins = str2double(partsStr(1)); if isnan(mins), mins=0; end
			secs = mins*60 + lastSec; return;
		elseif numel(partsStr)>=3
			hrs = str2double(partsStr(1)); mins = str2double(partsStr(2));
			if isnan(hrs), hrs=0; end; if isnan(mins), mins=0; end
			secs = hrs*3600 + mins*60 + lastSec; return;
		else
			secs = NaN; return;
		end
	end

	function doJumpTime(fig, timeStr)
		S = guidata(fig);
		if isempty(S) || isempty(S.tFull)
			errordlg('Please open an EDF file first.','No Data'); return;
		end
		t = S.tFull; fs = S.fs;
		sec = parseTimeStr(timeStr);
		if isempty(sec) || isnan(sec)
			beep; warndlg('Failed to parse time. Please enter seconds or hh:mm:ss format.','Invalid time format'); return;
		end
		% Clamp target time within bounds
		target = max(t(1), min(t(end), sec));
		halfW = max(1/fs, S.jumpWindowSec/2);
		x1 = max(t(1), target - halfW);
		x2 = min(t(end), target + halfW);
		if x2 <= x1, x2 = min(t(end), x1 + max(1/fs, S.jumpWindowSec)); end
		% Window indices
		i1 = max(1, floor(x1*fs)+1);
		i2 = min(numel(t), ceil(x2*fs)+1);
		% Compute Y limits and pad
		if ~isempty(S.ecgDispFull)
			ySegTop = S.ecgDispFull(i1:i2);
			yMinT = min(ySegTop); yMaxT = max(ySegTop);
			if yMinT==yMaxT, yMinT=yMinT-1; yMaxT=yMaxT+1; end
			padT = 0.1*(yMaxT - yMinT + eps);
			yLimTop = [yMinT - padT, yMaxT + padT];
		else
			yLimTop = get(S.axTop,'YLim');
		end
		if ~isempty(S.ecgRaw)
			ySegB = S.ecgRaw(i1:i2);
			yMinB = min(ySegB); yMaxB = max(ySegB);
			if yMinB==yMaxB, yMinB=yMinB-1; yMaxB=yMaxB+1; end
			padB = 0.1*(yMaxB - yMinB + eps);
			yLimB = [yMinB - padB, yMaxB + padB];
		else
			yLimB = get(S.axBottom,'YLim');
		end
		set(S.axTop, 'XLim',[x1 x2], 'YLim', yLimTop);
		set(S.axBottom, 'XLim',[x1 x2], 'YLim', yLimB);
		refreshView(fig);
	end

	function v = safeField(st, fn, def)
		if isstruct(st) && isfield(st, fn) && ~isempty(st.(fn))
			v = st.(fn);
		else
			v = def;
		end
	end

	function idx = sanitizeIndices(idx, maxLen)
		idx = round(double(idx(:)'));
		idx = idx(idx>=1 & idx<=maxLen);
	end

	function deleteValidHandles(h)
		if isempty(h), return; end
		for kk = 1:numel(h)
			if ishghandle(h(kk))
				try
					delete(h(kk));
				catch
				end
			end
		end
	end

	function hTxt = drawTypeLabels(ax, xVals, yVals, typeList, offsetFrac)
		if isempty(xVals)
			hTxt=gobjects(0);
			return;
		end
		xl = get(ax,'XLim'); yl=get(ax,'YLim'); dx=(xl(2)-xl(1))*offsetFrac(1); dy=(yl(2)-yl(1))*offsetFrac(2);
		n = numel(xVals); hTxt = gobjects(n,1);
		for ii=1:n
			lbl = typeList{ii}; if ~ischar(lbl), lbl = char(string(lbl)); end
			try
				hTxt(ii) = text(ax, xVals(ii)+dx, yVals(ii)+dy, lbl, 'Color',[0.85 0 0], ...
					'FontSize',8, 'Clipping','on', 'Interpreter','none', 'VerticalAlignment','bottom');
			catch
			end
		end
	end

	function [rIdx, rTypes, annFs, annPath, rSecs] = tryLoadAnnotation(edfPath)
		rIdx=[]; rTypes={}; annFs=[]; annPath=''; rSecs=[];
		[edfDir, edfBase, ~] = fileparts(edfPath);
		[annRoot, subRel] = inferAnnotationsRoot(edfDir);
		if isempty(annRoot) || ~isfolder(annRoot), return; end
		cands = {
			fullfile(annRoot, subRel, [edfBase '-rpoint.csv'])
			fullfile(annRoot, subRel, [edfBase '-rpoints.csv'])
			fullfile(annRoot, [edfBase '-rpoint.csv'])
			fullfile(annRoot, [edfBase '-rpoints.csv'])
		};
		csvPath=''; for ii=1:numel(cands), if exist(cands{ii},'file'), csvPath=cands{ii}; break; end; end
		if isempty(csvPath)
			d1=dir(fullfile(annRoot,'**',[edfBase '-rpoint.csv'])); if isempty(d1), d1=dir(fullfile(annRoot,'**',[edfBase '-rpoints.csv'])); end
			if ~isempty(d1), csvPath=fullfile(d1(1).folder, d1(1).name); end
		end
		if isempty(csvPath), return; end
		try
			T=readtable(csvPath);
		catch
			return;
		end
		vn = lower(T.Properties.VariableNames);
		ir = find(strcmp(vn,'rpoint'),1); if isempty(ir), ir = find(contains(vn,'rpoint'),1); end
		it = find(strcmp(vn,'type'),1);
		isec = find(strcmp(vn,'seconds'),1);
		isr = find(strcmp(vn,'samplingrate'),1);
		if ~isempty(isec)
			try
				rSecs = double(T{:,isec});
			catch
				rSecs = [];
			end
		end
		if ~isempty(it)
			if iscell(T{:,it}), rTypes = T{:,it}; else, rTypes = cellstr(string(T{:,it})); end
		end
		if isempty(rSecs) && ~isempty(ir)
			try
				rIdx = T{:,ir};
			catch
				rIdx = [];
			end
		end
		if ~isempty(isr)
			annFs = T{1,isr}; if istable(annFs)||iscell(annFs), annFs = str2double(string(annFs)); end
			annFs = double(annFs); if isnan(annFs)||annFs<=0, annFs=[]; end
		end
		if isempty(annFs), annFs = S.fs_default; end
		annPath = csvPath;
	end

	function [annRoot, subRel] = inferAnnotationsRoot(edfDir)
		subRel=''; parts=strsplit(edfDir, filesep);
		idxEdps = find(strcmpi(parts,'edfs'),1,'last'); idxPoly = find(strcmpi(parts,'polysomnography'),1,'last');
		if isempty(idxPoly), annRoot=''; return; end
		annBaseParts = parts(1:idxPoly); annRoot = fullfile(annBaseParts{:}, 'annotations-rpoints');
		if ~isempty(idxEdps) && idxEdps < numel(parts)
			subRel = fullfile(parts{idxEdps+1:end});
		end
	end

	function rIdxOut = mapIndicesToFs(rIdxIn, fsIn, fsOut, maxLen)
		rIdxOut = round((double(rIdxIn)-1) * (fsOut/fsIn)) + 1;
		rIdxOut(rIdxOut<1)=1; rIdxOut(rIdxOut>maxLen)=maxLen; rIdxOut=unique(rIdxOut(:)');
	end

end
