function runEEGAnalysis(ana)
ft_defaults;
ana.table.Data			=[]; 
ana.warning.Color		= [ 0.5 0.5 0.5 ];
ana.codeVersion			= '1.10';
ana.versionLabel.Text	= [ana.versionLabel.UserData ' Code: V' ana.codeVersion];
colours					= analysisCore.optimalColours(10);
info					= load(ana.MATFile);
info.origSettings		= info.ana; info = rmfield(info,'ana');
info.sM=[]; info.eT=[]; info.origSettings.cal=[];info.origSettings.outcal=[];
info.seq.getLabels();
vars					= getVariables();
data_raw = []; trl=[]; triggers=[]; events=[]; timelock = []; freq = [];
fprintf('\n\n\n=====>>>>> EEG Analysis Initiated... <<<<<=====\n'); ts=tic;

%=====================================================PLOT RAW DATA TRIGGERS
if ana.plotTriggers
	info.seq.showLog(); drawnow;
	cfgRaw				= [];
	cfgRaw.dataset		= ana.EDFFile;
	cfgRaw.header		= ft_read_header(cfgRaw.dataset); 
	if contains(cfgRaw.header.chanunit{1},'mV'); mV = 1e3; else; mV = 0; end
	disp('============= HEADER INFO, please check! ====================');
	disp(cfgRaw.header); disp(cfgRaw.header.orig);
	if cfgRaw.header.nChans ~= ana.pDiode
		ana.pDiode		= cfgRaw.header.nChans;
		ana.bitChannels = ana.pDiode-8:ana.pDiode-1;
		ana.dataChannels = 1:ana.bitChannels(1)-1;
		warndlg('GUI channel assignments are incorrect, will correct this time!')
	end
	cfgRaw.continuous	= 'yes';
	cfgRaw.channel		= 'all';
	cfgRaw.demean		= 'yes';
	cfgRaw.detrend		= 'yes';
	cfgRaw.polyremoval  = 'yes';
	cfgRaw.chanindx     = ana.bitChannels;
	cfgRaw.threshold	= ana.threshold;
	cfgRaw.jitter		= ana.jitter;
	cfgRaw.minTrigger	= ana.minTrigger;
	cfgRaw.preTime		= ana.preTime;
	cfgRaw.correctID	= ana.correctID;
	data_raw			= ft_preprocessing(cfgRaw);
	cfgRaw.denoise		= false;
	[trl, events, triggers] = loadCOGEEG(cfgRaw);
	if isempty(trl)
		fprintf('--->>> NO Trials loaded\n');
	else
		fprintf('--->>> %i Trials loaded, plotting...\n',size(trl,1));
	end
	plotRawChannels(); drawnow;
	if ~isempty(trl) && size(trl,2) == 4
		plotTable(info.seq.outIndex,trl(:,4));
	end
	info.data_raw		= data_raw;
	info.events			= events;
	info.triggers		= triggers;
	info.trl			= trl;
	assignin('base','info',info);
	return; %we don't do any further analysis
end

%====================================================PARSE DATA AS TRIALS
cfg						= [];
cfg.dataset				= ana.EDFFile;
cfg.header				= ft_read_header(cfg.dataset); disp(cfg.header);
if contains(cfg.header.chanunit{1},'mV'); mV = 1e3; else; mV = 1; end %convert to microvolts
if cfg.header.nChans ~= ana.pDiode
	ana.pDiode			= cfg.header.nChans;
	ana.bitChannels		= ana.pDiode-8:ana.pDiode-1;
	ana.dataChannels	= 1:ana.bitChannels(1)-1;
	warndlg('GUI channel assignments were incorrect, will correct this time!')
end
cfg.continuous			= 'yes';
cfg.trialfun			= 'loadCOGEEG';
cfg.chanindx			= ana.bitChannels;
cfg.threshold			= ana.threshold;
cfg.jitter				= ana.jitter;
cfg.minTrigger			= ana.minTrigger;
cfg.correctID			= ana.correctID;
cfg.preTime				= ana.preTime;
cfg.denoise				= false;
cfg						= ft_definetrial(cfg); %find trials
cfg						= rmfield(cfg,'denoise');
cfg.demean				= ana.demean;

if strcmpi(ana.demean,'yes') 
	cfg.baselinewindow	= ana.baseline;
end
cfg.medianfilter		= ana.medianfilter;
cfg.dftfilter			= ana.dftfilter;
cfg.dftfreq				= [50 100 150];
cfg.detrend				= ana.detrend;
cfg.polyremoval			= ana.polyremoval;
cfg.channel				= ana.dataChannels;
if ana.rereference > 0 && any(ana.dataChannels == ana.rereference)
	cfg.reref			= 'yes';
	cfg.refchannel		= cfg.header.label{ana.rereference};
	cfg.refmethod		= ana.rerefMethod;
end
if ana.bandpass == true && ana.lowpass > 0 && ana.highpass > 0
	cfg.bpfilter		= 'yes';
	cfg.bpfreq			= [ana.highpass ana.lowpass];
	cfg.bpfilttype		= ana.filtertype;
	cfg.bpfiltdir		= ana.direction;
	if ana.order > 0; cfg.bpfiltord = ana.order; end
else
	if ana.lowpass > 0 
		cfg.lpfilter	= 'yes';
		cfg.lpfreq		= ana.lowpass;
		cfg.lpfilttype	= ana.filtertype;
		cfg.lpfiltdir   = ana.direction;
		if ana.order > 0; cfg.lpfiltord = ana.order; end
	end
	if ana.highpass > 0
		cfg.hpfilter	= 'yes';
		cfg.hpfreq		= ana.highpass;
		cfg.hpfilttype	= ana.filtertype;
		cfg.hpfiltdir   = ana.direction;
		if ana.order > 0; cfg.hpfiltord = ana.order; end
	end
end
if ana.plotFilter; cfg.plotfiltresp = 'yes'; end
data_eeg				= ft_preprocessing(cfg); %Load and filter data

if ana.makeSurrogate
	ts=tic;
	makeSurrogate(); %create some artificial data
	fprintf('\n===>>> Surrogate data took %.2f secs to generate...\n',toc(ts));
end

info.rejected = [];
if ana.rejectvisual % visual reject and load/save selected trials to GUI
	cfg					= [];
	cfg.box				= 'yes';
	cfg.latency			= 'all';
	cfg.method			= ana.rejecttype;
	index				= 1:length(data_eeg.trial);
	if ~isempty(ana.rejecttrials)
		cfg.trials		= setdiff(index,ana.rejecttrials);
	end
	data_eeg			= ft_rejectvisual(cfg,data_eeg);
	if length(data_eeg.trial) < length(index)
		info.rejected	= setdiff(index,data_eeg.cfg.trials);
		disp(['Rejected Trials: ' num2str(info.rejected)]);
		ana.rjHandle.Value = regexprep(num2str(info.rejected),'\s+',' ');
	end
end

%================================================RUN TIMELOCK ANALYSIS
varmap					= unique(data_eeg.trialinfo);
timelock				= cell(length(varmap),1);
avgfn					= eval(['@nan' ana.avgmethod]);
if ana.doTimelock
	for jv = 1:length(varmap)
		cfg				= [];
		cfg.trials		= find(data_eeg.trialinfo==varmap(jv));
		cfg.covariance	= ana.tlcovariance;
		cfg.keeptrials	= ana.tlkeeptrials;
		cfg.removemean	= ana.tlremovemean;
		if ~isempty(ana.plotRange);	cfg.latency = ana.plotRange; end
		timelock{jv}	= ft_timelockanalysis(cfg,data_eeg);
	end
	
	plotTimeLock();
	plotFreqPower();
	plotCSF();
end

%================================================RUN TIMEFREQ ANALYSIS
freq					= cell(length(varmap),1);
if ana.doTimeFreq
	for jtf = 1:length(varmap)
		cfg				= [];
		cfg.trials		= find(data_eeg.trialinfo==varmap(jtf));
		cfg.channel		= 1;
		cfg.method		= 'mtmconvol';
		cfg.taper		= ana.freqtaper;
		cfg.pad			= 'nextpow2';
		cfg.foi			= ana.freqrange;                  % analysis 2 to 30 Hz in steps of 2 Hz
		cfg.t_ftimwin	= ones(length(cfg.foi),1).*0.2;   % length of time window = 0.5 sec
		if ~isempty(ana.plotRange) && isnumeric(ana.plotRange) && length(ana.plotRange)==2
			cfg.toi		= ana.plotRange(1):0.05:ana.plotRange(2);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
		else
			cfg.toi		= min(data_eeg.time{1}):0.05:max(data_eeg.time{1});
		end
		freq{jtf}		= ft_freqanalysis(cfg,data_eeg);
	end
	plotFrequency();
end

%info.timelock			= timelock;
%info.freq				= freq;
info.data_raw			= data_raw;
%info.data_eeg			= data_eeg;
info.triggers			= triggers;
info.analSettings		= ana;
assignin('base','info',info);
%assignin('base','timelock',timelock);

plotTable(info.seq.outIndex, data_eeg.trialinfo);
fprintf('===>>> Analysis took %.2f seconds\n', toc(ts));


%=============================================================================
%================================================================SUB FUNCTIONS
%=============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GET VARS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vars = getVariables()
	if isprop(info.seq,'varLabels') && ~isempty(info.seq.varLabels)
		vars			= info.seq.varLabels;
	else
		vars			= cell(info.seq.minBlocks,1);
		for i=1:length(vars)
			vars{i} = num2str(i);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTable(intrig,outtrig)
	col1 = intrig;if size(col1,1)<size(col1,2); col1=col1';end
	col2 = outtrig;if size(col2,1)<size(col2,2); col2=col2';end
	col3 = vars; if size(col3,1)<size(col3,2); col3=col3';end
	col4 = 1:length(col3); if size(col4,1)<size(col4,2); col4=col4';end
	
	if length(col1) ~= length(col2)
		warning('Input and output triggers are different!')
		ana.warning.Color = [ 0.8 0.3 0.3 ];
	else
		ana.warning.Color = [ 0.3 0.8 0.3 ];
	end

	maxn = max([length(col1) length(col2) length(col3) length(col4)]);
	if length(col1) < maxn; col1(end+1:maxn) = NaN; end
	if length(col2) < maxn; col2(end+1:maxn) = NaN; end
	if length(col3) < maxn
		col3 = [col3;repmat({''},maxn-length(col3),1)];
	end
	if length(col4) < maxn; col4(end+1:maxn) = NaN; end
	tdata = table(col1,col2,col3,col4,'VariableNames',{'Triggers Sent','Data Triggers','Stimulus Value','Index'});
	ana.table.Data = tdata;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT TIME LOCKED RESPONSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTimeLock()
	[p f e] = fileparts(ana.EDFFile);
	h = figure('Name',['TL Data: ' f e],'Units','normalized',...
		'Position',[0 0.025 0.25 0.9]);
	if length(timelock) > 8
		tl = tiledlayout(h,'flow','TileSpacing','compact');
	else
		tl = tiledlayout(h,length(timelock),1,'TileSpacing','compact');
	end
	mn = inf; mx = -inf;
	fprintf('\n--->>> Plotting Time-Locked Potentials: \n');
	for jj = 1:length(timelock)
		if jj == 1;fprintf(' #%03i...\n', jj);else; fprintf('\b\b\b\b\b\b\b\b\b #%03i...\n', jj);end
		nexttile(tl,jj); 
% 		if length(ana.tlChannels)>1
% 			cfg = [];
% 			cfg.linecolor = 'gbywrgbkywrgbkywrgbkywlrb';
% 			cfg.interactive = 'no';
% 			cfg.linewidth = 2;
% 			cfg.channel = ana.tlChannels;
% 			ft_singleplotER(cfg,timelock{jj});
% 		end
		hold on
		if isfield(timelock{jj},'avg')
			for i = 1:length(timelock{jj}.label)
				analysisCore.areabar(timelock{jj}.time, timelock{jj}.avg(i,:)*mV,...
					timelock{jj}.var(i,:)*mV,colours(i,:));
			end
		else
			for i = 1:length(timelock{jj}.label)
				h=plot(timelock{jj}.time',squeeze(timelock{jj}.trial(:,i,:))'*mV,...
					':','Color',colours(i,:),'DisplayName',timelock{jj}.label{i});
				for hl=1:length(h)
					h(hl).Annotation.LegendInformation.IconDisplayStyle = 'off';
					h(hl).UserData = ['Trial: ' num2str(timelock{jj}.cfg.trials(hl)) ' Ch: ' timelock{jj}.label{i}];
					h(hl).PickableParts = 'all';
					h(hl).ButtonDownFcn = @clickMe;
				end
				plot(timelock{jj}.time',avgfn(squeeze(timelock{jj}.trial(:,i,:)))'*mV,...
					'-','Color',colours(i,:),'LineWidth',1.5,'DisplayName',timelock{jj}.label{i});
			end
		end
		if isnumeric(ana.plotRange);xlim([ana.plotRange(1) ana.plotRange(2)]);end
		box on;grid on; grid minor; axis tight;
		if length(ana.tlChannels)>1 && jj == 1000
			legend(cat(1,{'AVG'},timelock{1}.label));
		elseif jj == 1
			legend(timelock{1}.label);
		end
		if min(ylim)<mn;mn=min(ylim);end
		if max(ylim)>mx;mx=max(ylim);end
		l = line([0 0],ylim,'LineStyle','--','LineWidth',1.25,'Color',[.4 .4 .4]);
		l.Annotation.LegendInformation.IconDisplayStyle = 'off';
		l.ButtonDownFcn = @cloneAxes;
		t = title(['Var: ' num2str(jj) ' = ' vars{jj}]);
		t.ButtonDownFcn = @cloneAxes;
		hz = zoom;hz.ActionPostCallback = @myCallbackZoom;
		hp = pan;hp.ActionPostCallback = @myCallbackZoom;
	end
	fprintf(' ... DONE\n');
	interv = info.origSettings.VEP.Flicker;
	nint = round(max(timelock{1}.time) / interv);
	for j = 1:length(timelock)
		nexttile(tl,j);
		ylim([mn mx]);
		for kk = 1:2:nint
			rectangle('Position',[(kk-1)*interv mn interv mx-mn],...
			'FaceColor',[0.8 0.8 0.8 0.1],'EdgeColor','none');
		end
	end
	t = sprintf('TL: dft=%s demean=%s (%.2f %.2f) detrend=%s poly=%s lp=%.2f hp=%.2f',ana.dftfilter,ana.demean,ana.baseline(1),ana.baseline(2),ana.detrend,ana.polyremoval,ana.lowpass,ana.highpass);
	tl.XLabel.String = 'Time (s)';
	tl.YLabel.String = 'Amplitude (\muv)';
	if ~isempty(info.rejected); t = [t '\newlineRejected Trials: ' num2str(info.rejected)]; end
	tl.Title.String = [t '\newlineComments: ' info.origSettings.comments];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT POWER ACROSS FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFreqPower()
	[~, f, e] = fileparts(ana.EDFFile);
	h = figure('Name',['TL Data: ' f e],'Units','normalized',...
		'Position',[0.25 0.025 0.25 0.9]);
	if length(timelock) > 8
		tl = tiledlayout(h,'flow','TileSpacing','compact');
	else
		tl = tiledlayout(h,length(timelock),1,'TileSpacing','compact');
	end
	ff = 1/info.origSettings.VEP.Flicker;
	mn = inf; mx = -inf;
	powf(length(timelock),1) = struct('f0',[],'f1',[],'f2',[],'A',[]);
	tlNames = {data_eeg.hdr.label{ana.tlChannels}};
	daT = table;
	daT.Properties.Description = ['TL Data: ' f e];
	daT{:,'Index'} = cell2mat(info.seq.varList(:,1));
	daT{:,'Labels'} = vars;
	daT{:,'Trials'} = info.seq.varList(:,2);
	for vi = 1 : info.seq.nVars
		daT{:,info.seq.nVar(vi).name} = cell2mat(info.seq.varList(:,vi+2));
	end
	fprintf('--->>> Plotting FFT Power for Condition: \n');
	for j = 1:length(timelock)
		powf(j).label = timelock{j}.label;
		powf(j).trials = timelock{j}.cfg.trials;
		if isfield(timelock,'trialinfo')
			powf(j).trialinfo = unique(timelock{j}.trialinfo);
		end
		if j == 1;fprintf(' #%03i...\n',j);else;fprintf('\b\b\b\b\b\b\b\b\b #%03i...\n',j);end
		minidx = analysisCore.findNearest(timelock{j}.time, ana.analRange(1));
		maxidx = analysisCore.findNearest(timelock{j}.time, ana.analRange(2));
		nexttile(tl,j)
		hold on
		for ch = 1:length(timelock{j}.label)
			if isfield(timelock{j},'avg')
				dt = timelock{j}.avg(ch,minidx:maxidx)*mV;
				[P,f,A,f0,f1,f2] = doFFT(dt);
				if ~ana.combineFFT
					plot(f,P,'Color',colours(ch,:));
				end
				if any(contains(tlNames,timelock{j}.label{ch}))
					powf(j).f0 = [powf(j).f0 f0];
					powf(j).f1 = [powf(j).f1 f1];
					powf(j).f2 = [powf(j).f2 f2];
					powf(j).A  = [powf(j).A rad2deg(A)];
				end
				if min(P)<mn;mn=min(P);end
				if max(P)>mx;mx=max(P);end
			else
				dt = squeeze(timelock{j}.trial(:,ch,minidx:maxidx))*mV;
				[P,f] = doFFT(avgfn(dt));
				if min(P)<mn;mn=min(P);end
				if max(P)>mx;mx=max(P);end
				h=plot(f,P,'--','Color',colours(ch,:));
				h.Annotation.LegendInformation.IconDisplayStyle = 'off';
				for jj = 1:size(dt,1)
					[P(jj,:),f,A,f0,f1,f2] = doFFT(dt(jj,:));
					if any(contains(tlNames,timelock{j}.label{ch}))
						powf(j).f0 = [powf(j).f0 f0];
						powf(j).f1 = [powf(j).f1 f1];
						powf(j).f2 = [powf(j).f2 f2];
						powf(j).A  = [powf(j).A rad2deg(A)];
					end
				end
				if ~ana.combineFFT
					[avg,err] = analysisCore.stderr(P, ana.errormethod, [], ana.pvalue, [], avgfn);
					analysisCore.areabar(f,avg,err,colours(ch,:));
				end
			end
			if ch == 1; PP = P; else; PP = [PP; P]; end
		end
		powf(j).PP = PP;
		if ana.combineFFT
			[avg,err] = analysisCore.stderr(PP, ana.errormethod, [], ana.pvalue, [], avgfn);
			analysisCore.areabar(f,avg,err,[1 0.5 0]);
		end
		daT{j,'Power'} = powf(j);
		l = line([[ff ff]',[ff*2 ff*2]'],[ylim' ylim'],'LineStyle','--','LineWidth',1.25,'Color',[.4 .4 .4]);
		l(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
		l(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
		if j==1 && ~ana.combineFFT; legend(timelock{1}.label);end
		box on;grid on; grid minor;
		t = title(['Var: ' num2str(j) ' = ' vars{j}]);
		t.ButtonDownFcn = @cloneAxes;
		hz = zoom;hz.ActionPostCallback = @myCallbackZoom;
		hp = pan;hp.ActionPostCallback = @myCallbackZoom;
	end
	fprintf(' ... DONE!\n');
	for jj = 1:length(timelock);nexttile(tl,jj);ylim([0 mx]);xlim([-1 35]);end
	t = sprintf('TL: dft=%s demean=%s (%.2f %.2f) detrend=%s poly=%s ANALTIME: %.2f-%.2f',ana.dftfilter,...
		ana.demean,ana.baseline(1),ana.baseline(2),ana.detrend,ana.polyremoval, ana.analRange(1), ana.analRange(2));
	tl.XLabel.String = 'Frequency (Hz)';
	if isfield(timelock{j},'avg')
		tl.YLabel.String = 'FFT Power \muv';
	else
		tl.YLabel.String = ['FFT Power \muv \pm' ana.errormethod];
	end
	if ~isempty(info.rejected); t = [t '\newlineRejected Trials: ' num2str(info.rejected)]; end
	tl.Title.String = [t '\newlineComments: ' info.origSettings.comments];
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TUNING CURVES
	[~, fL, eXT] = fileparts(ana.EDFFile);
	h = figure('Name',['F: ' fL eXT],'Units','normalized',...
		'Position',[0.2 0.2 0.6 0.6]);
	tl = tiledlayout(h,'flow','TileSpacing','compact');
	nexttile(tl)
	disp('--->>> Plotting Tuning Curves:'); warning off
	for i = 1:height(daT)
		if isfield(timelock{j},'avg')
			daT{i,'f0'} = avgfn(daT{i,'Power'}.f0);
			daT{i,'f0err'} = 0;
			daT{i,'f1'} = avgfn(daT{i,'Power'}.f1);
			daT{i,'f1err'} = 0;
			daT{i,'f2'} = avgfn(daT{i,'Power'}.f2);
			daT{i,'f2err'} = 0;
			daT{i,'A'} = avgfn(daT{i,'Power'}.A);
			daT{i,'Aerr'} = 0;
		else
			[daT{i,'f0'},e] = analysisCore.stderr(daT{i,'Power'}.f0, ana.errormethod, [], ana.pvalue, [], avgfn);
			daT{i,'f0err'} = e';
			[daT{i,'f1'},e] = analysisCore.stderr(daT{i,'Power'}.f1, ana.errormethod, [], ana.pvalue, [], avgfn);
			daT{i,'f1err'} = e';
			[daT{i,'f2'},e] = analysisCore.stderr(daT{i,'Power'}.f2, ana.errormethod, [], ana.pvalue, [], avgfn);
			daT{i,'f2err'} = e';
			[daT{i,'A'},e] = analysisCore.stderr(daT{i,'Power'}.A, ana.errormethod, [], ana.pvalue, [], avgfn);
			daT{i,'Aerr'}  = e';
		end
	end
	xa = 1:height(daT);
	if info.seq.addBlank
		xb = [ daT.Index(end); daT.Index(1:end-1) ];
		daT = daT(xb,:);
		daT.Order = xa';
		for jj = 1:height(daT)
			if jj == 1
				xlab{jj} = ['blank:' num2str(xb(jj))];
			else
				xlab{jj} = num2str(xb(jj));
			end
		end
	else
		daT.Order = xa';
		for jj = 1:length(xb)
			xlab{jj} = num2str(xb(jj));
		end
	end
	daT.Properties.RowNames = xlab;
	daT = movevars(daT,'Order','After','Index');
	
	if length(daT{1,'f1err'})==2
		thrsh = daT{1,'f1err'}(1,2);
	else
		if isempty(regexpi(ana.errormethod,'SE|SD'))
			thrsh = daT{1,'f1'} + daT{1,'f1err'};
		else
			thrsh = daT{1,'f1'} + daT{1,'f1err'}*2;
		end
	end
	info.daT = daT;
	opts = {'Marker','.','MarkerSize',16,'LineStyle','-'};
	if max(daT.f0err)==0
		yyaxis right
		plot(daT.Order,daT.A,'o--'); ylabel('Phase (deg)')
		yyaxis left
		pl = plot(daT.Order,[daT.f0 daT.f1 daT.f2],opts{:});
		pl(1).Color = colours(1,:);pl(2).Color = colours(2,:);pl(3).Color = colours(3,:);
		pl(1).Parent.XTick = daT.Order;
		pl(1).Parent.XTickLabel = daT.Properties.RowNames;
		pl(1).Parent.XTickLabelRotation=45;
	else
		hold on
		yyaxis right
		analysisCore.areabar(daT.Order,daT.A,daT.Aerr,[1 0.5 0],0.2); ylabel('Phase (deg)')
		yyaxis left
		pl = analysisCore.areabar(daT.Order,daT.f0,daT.f0err,colours(1,:),0.1,opts{:});
		pl = analysisCore.areabar(daT.Order,daT.f1,daT.f1err,colours(2,:),0.25,opts{:});
		pl = analysisCore.areabar(daT.Order,daT.f2,daT.f2err,colours(3,:),0.1,opts{:});
		pl = pl.plot;
		pl(1).Parent.XTick = daT.Order;
		pl(1).Parent.XTickLabel = daT.Properties.RowNames';
		pl(1).Parent.XTickLabelRotation=45;
		lp = line(xlim, [thrsh thrsh],'LineStyle','--','LineWidth',2,'Color',[.9 0 0]);
		lp.Annotation.LegendInformation.IconDisplayStyle = 'off';
	end
	xlim([0.9 length(xa)+0.1]);
	ymax = max(ylim);
	legend({'0th','1st','2nd'});box on; grid on;
	title(['Flicker Frequency: ' num2str(ff) 'Hz'])
	xlabel('Variable #');
	
	if info.seq.nVars == 2
		v1 = []; v2 = []; v1.threshold = thrsh; v2.threshold = thrsh;
		v1.name = daT.Properties.VariableNames{5};
		v2.name = daT.Properties.VariableNames{6};
		v1.values = unique(daT.(v1.name)); v1.values(isnan(v1.values)) = [];
		v2.values = unique(daT.(v2.name)); v2.values(isnan(v2.values)) = [];
		v1.n = length(v1.values);
		v2.n = length(v2.values);
		if contains(daT.Properties.RowNames{1},'blank')
			v1.hasBlank = true; v2.hasBlank = true;
		else
			v1.hasBlank = false; v2.hasBlank = false;
		end
		
		for jj = 1:v1.n
			if contains(daT.Properties.RowNames{1},'blank')
				v1.idx{jj} = [1];
				v1.label{jj} = [0];
			else
				v1.idx{jj} = [];
				v1.label{jj} = [];
			end
			v1.idx{jj} = [v1.idx{jj}; find(daT.(v1.name) == v1.values(jj))];
			v1.label{jj} = [v1.label{jj}; v2.values];
		end
		
		for jj = 1:v2.n
			if contains(daT.Properties.RowNames{1},'blank')
				v2.idx{jj} = [1];
				v2.label{jj} = [0];
			else
				v2.idx{jj} = [];
				v2.label{jj} = [];
			end
			v2.idx{jj} = [v2.idx{jj}; find(daT.(v2.name) == v2.values(jj))];
			v2.label{jj} = [v2.label{jj}; v1.values];
		end
		
		v1.x = 1:length(v1.idx{1});
		v1.nOther = v2.n;
		v1.nameOther = v2.name;
		if v1.hasBlank; v1.labels = [0 v2.values']; else; v1.labels = [v2.values']; end
		v2.x = 1:length(v2.idx{1});
		v2.nOther = v1.n;
		v2.nameOther = v1.name;
		if v2.hasBlank; v2.labels = [0 v1.values']; else; v2.labels = [v1.values']; end
		for xx = 1 : 2
			l=[];
			if xx == 1
				pI = v2; pO = v1; 
			else
				pI = v1; pO = v2; 
				h = figure('Name',['F: ' fL eXT],'Units','normalized',...
					'Position',[0.1 0.1 0.8 0.8]);
				tl = tiledlayout(h,'flow','TileSpacing','compact');
			end
			for jj = 1 : length(pI.idx)
				if jj == 1;fprintf(' #%03i...\n', jj);else; fprintf('\b\b\b\b\b\b\b\b\b #%03i...\n', jj);end
				nexttile(tl); hold on
				if ana.linearAxis == true
					x = pI.x;
				else
					x = pI.labels;
				end
				if max(max(daT.f0err))==0
					yyaxis right;
					plot(x,daT.A(pI.idx{jj})','o--');ylabel('Phase (deg)')
					yyaxis left;
					if xx == 1
						points=[daT.f0(pI.idx{jj})'; daT.f1(pI.idx{jj})'; daT.f2(pI.idx{jj})']';
					else
						points=daT.f1(pI.idx{jj});
					end
					pl = plot(x,points,opts{:});
					if xx == 1
						pl(1).Color = colours(1,:);pl(2).Color = colours(2,:);pl(3).Color = colours(3,:);
					else
						pl(1).Color = colours(2,:);
					end
					if contains(pI.nameOther,'contrast')
						if pI.hasBlank && ana.removeBlank
							y=daT.f1(pI.idx{jj}) - daT.f1(1);
						else
							y=daT.f1(pI.idx{jj});
						end
						if ~isempty(ana.excludePoints) && ~contains(ana.excludePoints,'none')
							midx = eval(ana.excludePoints);
						else
							midx = 1:length(x);
						end
						l{jj} = fitlm(x(midx),y(midx),'linear','VarNames',{'contrast','power'});
						plot(l{jj});
						anova(l{jj},'summary')
					end
				else
					yyaxis right
					analysisCore.areabar(x,daT.A(pI.idx{jj}),daT.Aerr(pI.idx{jj},:),[1 0.5 0],0.1);ylabel('Angle (deg)')
					yyaxis left
					if xx == 1
						pl = analysisCore.areabar(x, daT.f0(pI.idx{jj}), daT.f0err(pI.idx{jj},:),colours(1,:),0.07,opts{:});
						pl = analysisCore.areabar(x, daT.f2(pI.idx{jj}), daT.f2err(pI.idx{jj},:),colours(3,:),0.07,opts{:});
					end
					pl = analysisCore.areabar(x, daT.f1(pI.idx{jj}), daT.f1err(pI.idx{jj},:),colours(2,:),0.2,opts{:});
					pl = pl.plot;
					lp = line([x(1) x(end)], [thrsh thrsh],'LineStyle','--','LineWidth',2,'Color',[.9 0 0]);
					lp.Annotation.LegendInformation.IconDisplayStyle = 'off';
					if contains(pI.nameOther,'contrast')
						if pI.hasBlank && ana.removeBlank
							y=daT.f1(pI.idx{jj}) - daT.f1(1);
						else
							y=daT.f1(pI.idx{jj});
						end
						if ~isempty(ana.excludePoints) && ~contains(ana.excludePoints,'none')
							midx = eval(ana.excludePoints);
						else
							midx = 1:length(x);
						end
						l{jj} = fitlm(x(midx),y(midx),'linear','VarNames',{'contrast','power'});
						info.fitX{jj} = x(midx);
						info.fitY{jj} = y(midx);
						plot(l{jj});
						anova(l{jj},'summary')
					end
				end
				
				pl(1).Parent.XTick = x;
				pl(1).Parent.XTickLabel = pI.label{jj};
				pl(1).Parent.XTickLabelRotation=45;
				xlim([x(1) - (x(end)/20) x(end) + (x(end)/20)]);
				ylim([-inf ymax]);
				t=title(['Power at ' pI.name ': ' num2str(pI.values(jj))]);
				t.ButtonDownFcn = @cloneAxes;
				xlabel(pO.name);
				box on;grid on; grid minor;
			end
			if ~isempty(l);	info.fits = l; end
		end
	end
	info.v1 = v1;
	info.v2 = v2;
	if isfield(timelock{j},'avg')
		tl.YLabel.String = 'FFT Power';
	else
		tl.YLabel.String = ['FFT Power \pm' ana.errormethod];
	end
	tl.Title.String = ['Tuning Averages for Channels ' num2str(ana.tlChannels)];
	figure(h);drawnow
	fprintf(' ... DONE!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT TIME FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCSF()
	
	
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT TIME FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFrequency()
	[~, f, e] = fileparts(ana.EDFFile);
	h = figure('Name',['TF Data: ' f e],'Units','normalized',...
		'Position',[0.6 0.025 0.25 0.9]);
	if length(freq) > 8
		tl = tiledlayout(h,'flow','TileSpacing','compact');
	else
		tl = tiledlayout(h,length(freq),1,'TileSpacing','compact');
	end
	
	for jj = 1:length(freq)
		nexttile(tl);
		cfg = [];
		if ~contains(ana.freqbaseline,'none')
			cfg.baseline = ana.freqbaselinevalue;
			cfg.baselinetype = ana.freqbaseline;
		end
		ft_singleplotTFR(cfg,freq{jj});
		line([0 0],[min(ana.freqrange) max(ana.freqrange)],'LineWidth',2);
		xlabel('Time (s)');
		ylabel('Frequency (Hz)');
		box on;grid on; axis tight
		t =title(['Var: ' num2str(jj) ' = ' vars{jj}]);
		t.ButtonDownFcn = @cloneAxes;
	end
	tl.Title.String = 'Time Frequency Analysis';	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRawChannels()
	% plotting code to visualise the raw data triggers
	offset = 0;
	nchan = length(cfgRaw.header.label);
	h = figure('Name',['RAW Data: ' cfgRaw.dataset],'Units','normalized',...
		'Position',[0.05 0.05 0.4 0.9]);
	tl = tiledlayout(h,nchan,1,'TileSpacing','compact','Padding','none');
	tm = data_raw.time{1};
    if ~isempty(trl)
        xl = [tm(trl(1,1))-1 tm(trl(1,1))+9];
    else
        xl = [10 20];
    end
	for i = 1:nchan
		ch{i} = data_raw.trial{1}(i+offset,:);
		baseline = nanmedian(ch{i});
		ch{i} = (ch{i} - baseline);
		ch{i} = ch{i} / max(ch{i});
		nexttile(tl,i)
		p = plot(tm,ch{i},'k-');
		dtt = p.DataTipTemplate;
		dtt.DataTipRows(1).Format = '%.3f';
		line([min(tm) max(tm)], [0 0],'LineStyle',':','Color',[0.4 0.4 0.4]);
		hold on
		if ~any(ana.bitChannels == i) && (i == 1 || i == ana.pDiode)
			for ii = 1:length(events)
				if ~isempty(events(ii).times)
					y = repmat(ii/10, [1 length(events(ii).times)]);
					plot(events(ii).times,y,'.','MarkerSize',12);
				end
			end
			ylim([-inf inf]);
		elseif any(ana.bitChannels == i)
			ii = i - (ana.bitChannels(1)-1);
			if ~isempty(events(ii).times)
				p=plot(events(ii).times,0.75,'r.','MarkerSize',12);
				dtt = p.DataTipTemplate;
				dtt.DataTipRows(1).Format = '%.3f';
			end
			ylim([-0.05 1.05]);
		end
		if any([ana.dataChannels ana.pDiode] == i) && i == 1 && ~isempty(trl) && size(trl,1) > 1
			ypos = 0.2;
			for jj = 1:size(trl,1) 
				line([tm(trl(jj,1)) tm(trl(jj,2))],[ypos ypos]);
				plot([tm(trl(jj,1)) tm(trl(jj,1)-trl(jj,3)) tm(trl(jj,2)+trl(jj,3))],ypos,'ko','MarkerSize',8);
				%text(tm(trl(jj,1)-trl(jj,3)),ypos,['\leftarrow' num2str(trl(jj,4))]);
				%text(tm(trl(jj,2)+trl(jj,3)),ypos,'\leftarrow255');
				ypos = ypos+0.125;
				if ypos > 1.0; ypos = 0.3;end
			end
			trgVals = num2cell([triggers.value]);
			trgVals = cellfun(@num2str,trgVals,'UniformOutput',false);
			trgTime = [triggers.time];
			trgY = ones(1,length(trgTime));
			text(trgTime,trgY,trgVals);
		end
		title(data_raw.label{i});
		xlim(xl);
	end
	hz = zoom;
	hz.ActionPostCallback = @myCallbackScroll;
	hp = pan;
	hp.enable = 'on';
	hp.Motion = 'horizontal';
	hp.ActionPostCallback = @myCallbackScroll;
	tl.XLabel.String = 'Time (s)';
	tl.YLabel.String = 'Normalised Amplitude';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCallbackScroll(~,event)
	src = event.Axes;
	xl = src.XLim;
	for i = 1:length(src.Parent.Children)
		if i < length(src.Parent.Children)-2
			src.Parent.Children(i).YLim = [-0.05 1.05];
		else
			ylim(src.Parent.Children(i),'auto');
		end
		if ~all(xl == src.Parent.Children(i).XLim)
			src.Parent.Children(i).XLim = xl;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCallbackZoom(~,event)
	src = event.Axes;
	xl = src.XLim;
	xy = src.YLim;
	for i = 1:length(src.Parent.Children)
		if isa(src.Parent.Children(i),'matlab.graphics.axis.Axes')
			if ~all(xl == src.Parent.Children(i).XLim)
				src.Parent.Children(i).XLim = xl;
			end
			if ~all(xy == src.Parent.Children(i).YLim)
				src.Parent.Children(i).YLim = xy;
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clickMe(src, ~)
	if ~exist('src','var')
		return
	end
	
	if src.LineWidth < 1
		src.LineWidth = 1; src.LineStyle =  '-';
	else
		src.LineWidth = 0.5; src.LineStyle =  ':';
	end
	ud = get(src,'UserData');
	if ~isempty(ud)
		disp(ud);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cloneAxes(src,~)
	disp('Cloning axis!')
	if ~isa(src,'matlab.graphics.axis.Axes')
		if isa(src.Parent,'matlab.graphics.axis.Axes')
			src = src.Parent;
		end
	end
	f=figure;
	nsrc = copyobj(src,f);
	nsrc.OuterPosition = [0.05 0.05 0.9 0.9];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, f, A, p0, p1, p2] = doFFT(p,fs,ff,normalise,useHanning)
	if ~exist('fs','var')||isempty(fs); fs = data_eeg.fsample; end
	if ~exist('ff','var')||isempty(ff); ff = (1/info.origSettings.VEP.Flicker); end
	if ~exist('normalise','var')||isempty(normalise); normalise = true; end
	if ~exist('useHanning','var')||isempty(useHanning); useHanning = ana.fftwindow; end
	
	L = length(p);
	
	if useHanning
		win = hanning(L, 'periodic');
		Pi = fft(p.*win'); 
	else
		Pi = fft(p);
	end

	if normalise
		P2 = abs(Pi/L);
		P=P2(1:floor(L/2)+1);
		P(2:end-1) = P(2:end-1)*2;
		f = fs * (0:(L/2))/L;
	else
		NumUniquePts = ceil((L+1)/2);
		P = abs(Pi(1:NumUniquePts));
		f = (0:NumUniquePts-1)*fs/L;
	end

	idx = analysisCore.findNearest(f, ff);
	p1 = P(idx);
	A = angle(Pi(idx));
	idx = analysisCore.findNearest(f, 0);
	p0 = P(idx);
	idx = analysisCore.findNearest(f, ff*2);
	p2 = P(idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doCSF()
	

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeSurrogate()
	randPhaseRange			= 0; %how much to randomise phase?
	rphase					= 0; %default phase
	basef					= 1; % base frequency
	onsetf					= 10; %an onset at 0 frequency
	onsetLength				= 3; %length of onset signal
	onsetDivisor			= 0.5; %scale the onset frequency
	burstf					= 30; %small burst frequency
	burstOnset				= 1.0; %time of onset of burst freq
	burstLength				= 1.0; %length of burst
	powerDivisor			= 2.0; %how much to attenuate the secondary frequencies
	group2Divisor			= 1; %do we use a diff divisor for group 2?
	noiseDivisor			= 1; %scale noise to signal
	piMult					= basef * 2; %resultant pi multiplier
	burstMult				= burstf * 2; %resultant pi multiplier
	onsetMult				= onsetf * 2; %onset multiplier
	lowpassNoise			= true;
	options = {['t|' num2str(randPhaseRange)], 'Random phase range in radians?';...
		['t|' num2str(rphase)], 'Default phase?';...
		['t|' num2str(basef)], 'Base Frequency (Hz)';...
		['t|' num2str(onsetf)], 'Onset (time=0) Frequency (Hz)';...
		['t|' num2str(onsetDivisor)], 'Onset F Power Divisor';...
		['t|' num2str(burstf)], 'Burst Frequency (Hz)';...
		['t|' num2str(burstOnset)], 'Burst Onset Time (s)';...
		['t|' num2str(burstLength)], 'Burst Length (s)';...
		['t|' num2str(powerDivisor)], 'Burst Power Divisor';...
		['t|' num2str(group2Divisor)], 'Burst Power Divisor for Group 2';...
		['t|' num2str(noiseDivisor)], 'Noise Divisor';...
		'x|¤Lowpass?','Filter noise?';...
		};
	answer = menuN('Select Surrogate options:',options);
	drawnow;
	if iscell(answer) && ~isempty(answer)
		randPhaseRange = eval(answer{1});
		rphase = str2num(answer{2});
		basef = str2num(answer{3});
		onsetf = str2num(answer{4});
		onsetDivisor = str2num(answer{5});
		burstf = str2num(answer{6});
		burstOnset = str2num(answer{7});
		burstLength = str2num(answer{8});
		powerDivisor = str2num(answer{9});
		group2Divisor = str2num(answer{10});
		noiseDivisor = str2num(answer{11});
		lowpassNoise = logical(answer{12});
	end

	f = data_eeg.fsample; 
	maxTime = max(data_eeg.time{1});
	if onsetLength > maxTime; onsetLength = maxTime - 0.1; end
	if burstLength > maxTime; burstLength = maxTime - 0.1; end

	for k = 1:length(data_eeg.trial)
		time = data_eeg.time{k};
		tLength = length(data_eeg.time{k});
		tmult = (tLength-1) / f; 
		mx = max(data_eeg.trial{k}(ana.surrogateChannel,:));
		mn = min(data_eeg.trial{k}(ana.surrogateChannel,:));
		rn = mx - mn;
		y = createSurrogate();
		y = y * rn; % scale to the voltage range of the original trial
		y = y + mn;
		data_eeg.trial{k}(ana.surrogateChannel,:) = y;
	end
	
	function y = createSurrogate()
		rphase = rand * randPhaseRange;
		%base frequency
		y = sin((0 : (pi*piMult)/f : (pi*piMult) * tmult)+rphase)';
		y = y(1:tLength);
		%burst frequency with different power in group 2 if present
		rphase = rand * randPhaseRange;
		yy = sin((0 : (pi*burstMult)/f : (pi*burstMult) * burstLength)+rphase)';
		if false
			yy = yy ./ group2Divisor;
		else
			yy = yy ./ powerDivisor;
		end
		%intermediate onset frequency
		rphase = rand * randPhaseRange;
		yyy = sin((0 : (pi*onsetMult)/f : (pi*onsetMult) * onsetLength)+rphase)';
		yyy = yyy ./ onsetDivisor;
		%find our times to inject yy burst frequency
		st = analysisCore.findNearest(time,burstOnset);
		en = st + length(yy)-1;
		if en > length(y); en = length(y); end
		y(st:en) = y(st:en) + yy(1:length(st:en));
		%add our fixed 0.4s intermediate onset freq
		st = analysisCore.findNearest(time,0);
		en = st + length(yyy)-1;
		if en > length(y); en = length(y); end
		y(st:en) = y(st:en) + yyy;
		%add our noise
		if lowpassNoise
			y = y + ((lowpass(rand(size(y)),300,f)-0.5)./noiseDivisor);
		else
			y = y + ((rand(size(y))-0.5)./noiseDivisor);
		end
		%normalise our surrogate to be 0-1 range
		y = y - min(y); y = y / max(y); % 0 - 1 range;
		%make sure we are a column vector
		if size(y,2) < size(y,1); y = y'; end
	end
end

end
