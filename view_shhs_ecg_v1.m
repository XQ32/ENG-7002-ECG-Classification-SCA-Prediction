function view_shhs_ecg_v1()
% File: view_shhs_ecg_v1.m
% Type: Function (GUI tool)
% Description:
%   Open an SHHS EDF ECG viewer. Display raw/filtered signals in linked axes and overlay rpoints annotations.
% Usage:
%   Run in MATLAB: view_shhs_ecg_v1
% Dependencies:
%   edfread, ecgFilter, and GUI callbacks/helpers within this file
% Maintainer: N/A  |  Version: 1.0  |  Date: 2025-08-26


	% GUI handles and state
	S = struct();
	S.fs_default = 125; % SHHS1 default sampling rate
	S.ecg = [];
	S.ecgFiltered = [];
	S.ecgFull = [];
	S.ecgFilteredFull = [];
	S.fs = S.fs_default;
	S.t = [];
	S.tFull = [];
	S.rIdx = [];
	S.rTypes = [];
	S.rIdxFull = [];
	S.rTypesFull = {};
	S.rSecsFull = [];
	S.annFs = [];
	S.edfPath = '';
	S.annPath = '';
	S.displayMaxPoints = 20000;    % max points per curve in one view
	S.displayMaxMarkers = 5000;    % max R-peak markers in one view
	S.maxLabels = 300;             % max text labels in one view
	S.labelOffsetFrac = [0.02 0.01];

	% Create UI
	S.fig = figure('Name','SHHS ECG Viewer', 'NumberTitle','off', 'Color','w', ...
		'Units','normalized', 'Position',[0.1 0.1 0.8 0.8]);

	% Menus (using classic figure Label/Callback properties)
	mFile = uimenu(S.fig, 'Label','File');
	uimenu(mFile, 'Label','Open EDF...', 'Accelerator','O', 'Callback', @(src,evt)onOpenEdf(src));
	uimenu(mFile, 'Label','Clear', 'Callback', @(src,evt)onClear(src));
	uimenu(mFile, 'Label','Exit', 'Callback', @(~,~)close(S.fig));

	% Two-axes layout
	S.tiled = tiledlayout(S.fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');
	S.axRaw = nexttile(S.tiled, 1);
	hold(S.axRaw, 'on'); grid(S.axRaw, 'on'); box(S.axRaw, 'on');
	title(S.axRaw, 'Raw ECG'); xlabel(S.axRaw, 'Time (s)'); ylabel(S.axRaw, 'Amplitude');
	S.axFilt = nexttile(S.tiled, 2);
	hold(S.axFilt, 'on'); grid(S.axFilt, 'on'); box(S.axFilt, 'on');
	title(S.axFilt, 'Raw ECG (Inverted)'); xlabel(S.axFilt, 'Time (s)'); ylabel(S.axFilt, 'Amplitude');

	% Plot handle placeholders
	S.hRaw = plot(S.axRaw, NaN, NaN, 'k-');
	S.hFilt = plot(S.axFilt, NaN, NaN, 'b-');
	S.hRRaw = plot(S.axRaw, NaN, NaN, 'rv', 'MarkerSize', 6, 'LineWidth', 1.0); % R-peak markers (raw)
	S.hRFilt = plot(S.axFilt, NaN, NaN, 'rv', 'MarkerSize', 6, 'LineWidth', 1.0); % R-peak markers (filtered)
	S.hTxtRaw = gobjects(0);
	S.hTxtFilt = gobjects(0);

	% Linked axes and zoom/pan callbacks
	linkaxes([S.axRaw, S.axFilt], 'x'); % Link X first, then sync Y in callback
	S.zoomObj = zoom(S.fig); S.zoomObj.Motion = 'both'; S.zoomObj.Enable = 'on';
	S.panObj = pan(S.fig); S.panObj.Enable = 'on';
	S.zoomObj.ActionPostCallback = @(obj,evd)onZoomPan(evd);
	S.panObj.ActionPostCallback = @(obj,evd)onZoomPan(evd);

	% Save state
	guidata(S.fig, S);

	% Internal: open EDF
	function onOpenEdf(src)
		fig = ancestor(src, 'figure');
		S = guidata(fig);
		startDir = fullfile(pwd, 'shhs', 'polysomnography', 'edfs');
		if ~isfolder(startDir)
			startDir = pwd;
		end
		[fileName, fileDir] = uigetfile({'*.edf','EDF Files (*.edf)'}, 'Select EDF file', startDir);
		if isequal(fileName,0)
			return;
		end
		edfPath = fullfile(fileDir, fileName);
		try
			TT = edfread(edfPath);
		catch ME
			errordlg(sprintf('edfread failed:\n%s', ME.message), 'Read Error');
			return;
		end

		% Find ECG channel
		varNames = TT.Properties.VariableNames;
		iEcg = find(strcmp(varNames,'ECG'), 1);
		if isempty(iEcg)
			lowNames = lower(varNames);
			iEcg = find(contains(lowNames,'ecg') | contains(lowNames,'ekg'), 1, 'first');
		end
		if isempty(iEcg)
			errordlg('ECG channel not found.', 'Missing Channel');
			return;
		end
		ecgVar = varNames{iEcg};

		% Flatten to column vector
		ecgCol = TT.(ecgVar);
		if iscell(ecgCol)
			try
				ecg = vertcat(ecgCol{:});
			catch
				ecg = [];
				for kk = 1:numel(ecgCol)
					v = ecgCol{kk};
					if isstring(v) || ischar(v)
						v = str2num(v); %#ok<ST2NM>
					end
					ecg = [ecg; v(:)]; %#ok<AGROW>
				end
			end
		elseif isnumeric(ecgCol)
			ecg = ecgCol(:);
		else
			errordlg(sprintf('ECG channel type not supported: %s', class(ecgCol)), 'Unsupported Type');
			return;
		end
		ecg = double(ecg);

		% Sampling rate: use default SHHS1=125; if annotations exist, align to their fs if needed
		fs = S.fs_default;

		% Filtering (SHHS1 already notched at 60 Hz -> power_line_freq=0; method 2)
		try
			ecgFiltered = ecgFilter(ecg, fs, 2, 0);
			% Viewer: invert filtered signal for display only
			ecgFiltered = -ecgFiltered;
		catch ME
			errordlg(sprintf('Filtering failed:\n%s', ME.message), 'Filtering Error');
			return;
		end

		% Time axis (seconds)
		t = (0:numel(ecg)-1)'/fs;

		% Try to load annotations (prefer seconds)
		[rIdx, rTypes, annFs, annPath, rSecs] = tryLoadAnnotation(edfPath);
		if ~isempty(rSecs)
			% Use seconds directly as time
			rTimes = rSecs(:)';
			% Corresponding indices (for nearest Y sampling)
			rIdxEcg = max(1, min(numel(ecg), round(rTimes*fs)+1));
		elseif ~isempty(rIdx)
			% Fallback: use indices + sampling-rate mapping
			rIdxEcg = mapIndicesToFs(rIdx, annFs, fs, numel(ecg));
			rTimes = (rIdxEcg-1)/fs;
		else
			rIdxEcg = [];
			rTimes = [];
			rTypes = {};
			annFs = [];
			annPath = '';
		end

		% Save full data
		S.ecgFull = ecg; S.ecgFilteredFull = ecgFiltered; S.fs = fs; S.tFull = t;
		S.ecg = []; S.ecgFiltered = []; S.t = [];
		S.rIdxFull = rIdxEcg; S.rSecsFull = rTimes; S.rTypesFull = rTypes; S.annFs = annFs; S.edfPath = edfPath; S.annPath = annPath;
		guidata(S.fig, S);

		% Initially show a 30 s window, then decimate according to view
		xStart = t(1); xEnd = min(t(end), xStart + 30);
		xlim(S.axRaw, [xStart xEnd]); xlim(S.axFilt, [xStart xEnd]);
		refreshView(S.fig);
	end

	% Internal: clear
	function onClear(src)
		fig = ancestor(src, 'figure');
		S = guidata(fig);
		set(S.hRaw, 'XData', NaN, 'YData', NaN);
		set(S.hFilt, 'XData', NaN, 'YData', NaN);
		set(S.hRRaw, 'XData', NaN, 'YData', NaN);
		set(S.hRFilt, 'XData', NaN, 'YData', NaN);
		deleteValidHandles(S.hTxtRaw); deleteValidHandles(S.hTxtFilt);
		S.hTxtRaw = gobjects(0); S.hTxtFilt = gobjects(0);
		S.ecg = []; S.ecgFiltered = []; S.t = []; S.fs = S.fs_default;
		S.ecgFull = []; S.ecgFilteredFull = []; S.tFull = [];
		S.rIdx = []; S.rTypes = []; S.rIdxFull = []; S.rTypesFull = {};
		S.annFs = []; S.edfPath = ''; S.annPath = '';
		guidata(S.fig, S);
	end

	% Internal: after zoom/pan, sync Y-axis and keep labels visible
	function onZoomPan(evd)
		fig = ancestor(evd.Axes, 'figure');
		S = guidata(fig);
		if isempty(S) || ~isfield(S,'axRaw') || ~ishandle(S.axRaw)
			return;
		end
		axSrc = evd.Axes;
		if axSrc == S.axRaw
			axDst = S.axFilt;
		else
			axDst = S.axRaw;
		end
		% Refresh view and sync YLimits (X is already linked by linkaxes)
		refreshView(fig);
		try
			set(axDst, 'YLim', get(axSrc,'YLim'));
		catch
			% ignore
		end
	end

	% Utility: draw type text labels
	function hTxt = drawTypeLabels(ax, xVals, yVals, typeList, offsetFrac)
		if isempty(xVals)
			hTxt = gobjects(0);
			return;
		end
		% Compute offsets relative to current axis ranges
		xl = get(ax,'XLim'); yl = get(ax,'YLim');
		dx = (xl(2)-xl(1))*offsetFrac(1);
		dy = (yl(2)-yl(1))*offsetFrac(2);
		n = numel(xVals);
		hTxt = gobjects(n,1);
		for ii = 1:n
			lbl = typeList{ii};
			if ~ischar(lbl)
				lbl = char(string(lbl));
			end
			try
				hTxt(ii) = text(ax, xVals(ii)+dx, yVals(ii)+dy, lbl, 'Color',[0.85 0 0], ...
					'FontSize',8, 'Clipping','on', 'Interpreter','none', 'VerticalAlignment','bottom');
			catch
				% If it fails, skip this label
			end
		end
	end

	% Utility: delete valid handles
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

	% Utility: map indices according to sampling rates
	function rIdxOut = mapIndicesToFs(rIdxIn, fsIn, fsOut, maxLen)
		rIdxOut = round((double(rIdxIn)-1) * (fsOut/fsIn)) + 1;
		rIdxOut(rIdxOut < 1) = 1;
		rIdxOut(rIdxOut > maxLen) = maxLen;
		rIdxOut = unique(rIdxOut(:)');
	end

	% Utility: try to locate and read annotation CSV (return raw indices/types and sampling rate)
	function [rIdx, rTypes, annFs, annPath, rSecs] = tryLoadAnnotation(edfPath)
		rIdx = []; rTypes = {}; annFs = []; annPath = ''; rSecs = [];
		[edfDir, edfBase, ~] = fileparts(edfPath);
		% Infer annotations-rpoints root directory
		[annRoot, subRel] = inferAnnotationsRoot(edfDir);
		if isempty(annRoot) || ~isfolder(annRoot)
			return;
		end
		% Candidate paths (prefer the same subdirectory structure as EDF)
		cands = {
			fullfile(annRoot, subRel, [edfBase '-rpoint.csv'])
			fullfile(annRoot, subRel, [edfBase '-rpoints.csv'])
			fullfile(annRoot, [edfBase '-rpoint.csv'])
			fullfile(annRoot, [edfBase '-rpoints.csv'])
		};
		csvPath = '';
		for ii = 1:numel(cands)
			if exist(cands{ii}, 'file')
				csvPath = cands{ii};
				break;
			end
		end
		% If not found, search recursively
		if isempty(csvPath)
			d1 = dir(fullfile(annRoot, '**', [edfBase '-rpoint.csv']));
			if isempty(d1)
				d1 = dir(fullfile(annRoot, '**', [edfBase '-rpoints.csv']));
			end
			if ~isempty(d1)
				csvPath = fullfile(d1(1).folder, d1(1).name);
			end
		end
		if isempty(csvPath)
			return;
		end
		try
			T = readtable(csvPath);
		catch
			return;
		end
		vn = lower(T.Properties.VariableNames);
		ir = find(strcmp(vn,'rpoint'),1);
		if isempty(ir)
			% Compatible with other header names
			candR = find(contains(vn,'rpoint'), 1);
			ir = candR;
		end
		it = find(strcmp(vn,'type'),1);
		isec = find(strcmp(vn,'seconds'),1);
		isr = find(strcmp(vn,'samplingrate'),1);
		% Prefer 'seconds' field
		if ~isempty(isec)
			try
				rSecs = double(T{:, isec});
			catch
				rSecs = [];
			end
		end
		% Types
		if ~isempty(it)
			if iscell(T{:, it})
				rTypes = T{:, it};
			else
				rTypes = cellstr(string(T{:, it}));
			end
		end
		% If there is no 'seconds', read 'rpoint' indices
		if isempty(rSecs) && ~isempty(ir)
			try
				rIdx = T{:, ir};
			catch
				rIdx = [];
			end
		end
		if ~isempty(isr)
			annFs = T{1, isr};
			if istable(annFs) || iscell(annFs), annFs = str2double(string(annFs)); end
			annFs = double(annFs);
			if isnan(annFs) || annFs <= 0
				annFs = [];
			end
		end
		if isempty(annFs), annFs = S.fs_default; end
		annPath = csvPath;
	end

	% Utility: infer annotations-rpoints root directory and subpath from EDF directory
	function [annRoot, subRel] = inferAnnotationsRoot(edfDir)
		subRel = '';
		parts = strsplit(edfDir, filesep);
		idxEdps = find(strcmpi(parts, 'edfs'), 1, 'last');
		idxPoly = find(strcmpi(parts, 'polysomnography'), 1, 'last');
		if isempty(idxPoly)
			annRoot = '';
			return;
		end
		if isempty(idxEdps) || idxEdps < idxPoly
			annBaseParts = parts(1:idxPoly);
			annRoot = fullfile(annBaseParts{:}, 'annotations-rpoints');
			return;
		end
		annBaseParts = parts(1:idxPoly);
		annRoot = fullfile(annBaseParts{:}, 'annotations-rpoints');
		if idxEdps < numel(parts)
			subRel = fullfile(parts{idxEdps+1:end});
		else
			subRel = '';
		end
	end

	% Refresh plot based on current view (decimation, limit object counts)
	function refreshView(fig)
		S = guidata(fig);
		if isempty(S) || isempty(S.tFull)
			return;
		end
		t = S.tFull; fs = S.fs;
		n = numel(t);
		xl = get(S.axRaw, 'XLim');
		x1 = max(t(1), xl(1));
		x2 = min(t(end), xl(2));
		if x2 <= x1
			x2 = min(t(end), x1 + 1/fs);
		end
		iStart = max(1, floor(x1*fs)+1);
		iEnd = min(n, ceil(x2*fs)+1);
		if iEnd <= iStart, iEnd = min(n, iStart+1); end

		% Decimate main curves
		segLen = iEnd - iStart + 1;
		step = max(1, ceil(segLen / S.displayMaxPoints));
		idx = iStart:step:iEnd;
		xp = t(idx);
		ypRaw = S.ecgFull(idx);
		ypFilt = S.ecgFilteredFull(idx);
		set(S.hRaw, 'XData', xp, 'YData', ypRaw);
		set(S.hFilt, 'XData', xp, 'YData', ypFilt);

		% R-peaks and labels (within current view only)
		if ~isempty(S.rIdxFull)
			mask = S.rIdxFull >= iStart & S.rIdxFull <= iEnd;
			rIdxV = S.rIdxFull(mask);
			rTypesV = S.rTypesFull(mask);
			rTimesV = S.rSecsFull(mask); % if seconds exist, use time
			if ~iscell(rTypesV)
				rTypesV = cellstr(string(rTypesV));
			end
			if ~isempty(rIdxV)
				mstep = max(1, ceil(numel(rIdxV) / S.displayMaxMarkers));
				rIdxV = rIdxV(1:mstep:end);
				rTypesV = rTypesV(1:mstep:end);
				if ~isempty(S.rSecsFull)
					rTimes = rTimesV(:)';
				else
					rTimes = (rIdxV-1)'/fs;
				end
				yRraw = S.ecgFull(rIdxV);
				yRfilt = S.ecgFilteredFull(rIdxV);
				set(S.hRRaw, 'XData', rTimes, 'YData', yRraw);
				set(S.hRFilt, 'XData', rTimes, 'YData', yRfilt);
				% Labels: apply count threshold
				deleteValidHandles(S.hTxtRaw); deleteValidHandles(S.hTxtFilt);
				if numel(rIdxV) <= S.maxLabels
					S.hTxtRaw = drawTypeLabels(S.axRaw, rTimes, yRraw, rTypesV, S.labelOffsetFrac);
					S.hTxtFilt = drawTypeLabels(S.axFilt, rTimes, yRfilt, rTypesV, S.labelOffsetFrac);
				else
					S.hTxtRaw = gobjects(0);
					S.hTxtFilt = gobjects(0);
				end
			else
				set(S.hRRaw, 'XData', NaN, 'YData', NaN);
				set(S.hRFilt, 'XData', NaN, 'YData', NaN);
				deleteValidHandles(S.hTxtRaw); deleteValidHandles(S.hTxtFilt);
				S.hTxtRaw = gobjects(0); S.hTxtFilt = gobjects(0);
			end
		else
			set(S.hRRaw, 'XData', NaN, 'YData', NaN);
			set(S.hRFilt, 'XData', NaN, 'YData', NaN);
			deleteValidHandles(S.hTxtRaw); deleteValidHandles(S.hTxtFilt);
			S.hTxtRaw = gobjects(0); S.hTxtFilt = gobjects(0);
		end

		drawnow limitrate
		guidata(fig, S);
	end
end


