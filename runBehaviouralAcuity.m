function runBehaviouralAcuity(ana)

global rM
if ~exist('rM','var') || isempty(rM)
	rM = arduinoManager();
end
rM.openGUI = false;
if ~rM.isOpen; open(rM); drawnow; end %open our reward manager


%===============compatibility for windows===================
%if ispc; PsychJavaTrouble(); end
KbName('UnifyKeyNames');

%===================Initiate out metadata===================
ana.date		= datestr(datetime);
ana.version		= Screen('Version');
ana.computer	= Screen('Computer');
ana.gpu			= opengl('data');
staircase = [];

%===================experiment parameters===================
ana.screenID	= max(Screen('Screens'));%-1;

%==========================================================================
%==================================================Make a name for this run
cd(ana.ResultDir)
ana.timeExp		= fix(clock());
if ~isempty(ana.subject)
	if ana.useStaircase; type = 'BASTAIR'; else type = 'BAMOC'; end %#ok<*UNRCH>
	nameExp = [type '_' ana.subject];
	c = sprintf(' %i',ana.timeExp);
	nameExp = [nameExp c];
	ana.nameExp = regexprep(nameExp,' ','_');
else
	ana.nameExp = 'debug';
end
cla(ana.plotAxis1);
cla(ana.plotAxis2);
cla(ana.plotAxis3);
cla(ana.plotAxis4);

nBlocks = ana.nBlocks;
nBlocksOverall = nBlocks * length(ana.contrastRange);

%==========================================================================
%===========================================================response values
YESBLANK = 1; YESTARGET = 2; UNSURE = 4; BREAKINIT = -100; BREAKBLANK = -10; 
BREAKTARGET = -1; BREAKEXCL = -5; UNDEFINED = 0;
saveMetaData();

%==========================================================================
%======================================================stimulus objects
% ---- blank disc.
blank = discStimulus();
blank.name = ['DISC' ana.nameExp];
blank.colour = [0.5 0.5 0.5];
blank.size = ana.discSize;
blank.sigma = ana.sigma;
% ---- target stimulus
target = discStimulus();
target.name = ['DISC' ana.nameExp];
target.colour = [0.5 0.5 0.5];
target.size = 1.5;
target.sigma = 5;
% ---- grat stimulus
grat = gratingStimulus();
grat.mask = true;
grat.useAlpha = true;
grat.sf = ana.SF;
grat.name = ['GRAT' ana.nameExp];
grat.size = blank.size;
grat.sigma = ana.sigma;
% ---- fixation cross
fixX = fixationCrossStimulus();
fixX.colour = [1 1 1];
fixX.colour2 = [0 0 0];
fixX.size = ana.spotSize;
fixX.alpha = ana.spotAlpha;
fixX.alpha2 = ana.spotAlpha;
fixX.lineWidth = ana.spotLine;
%----------combine them into a single meta stimulus
stimuli = metaStimulus();
stimuli.name = ana.nameExp;
stimuli{1} = blank;
stimuli{2} = target;
stimuli{3} = grat;
stimuli{4} = fixX;

%==========================================================================
%======================================================open the PTB screens
PsychDefaultSetup(2);
Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SkipSyncTests', 0);
Screen('Preference', 'DefaultFontSize',30);
sM						= screenManager();
sM.screen				= ana.screenID;
sM.windowed				= ana.windowed;
sM.pixelsPerCm			= ana.pixelsPerCm;
sM.distance				= ana.distance;
sM.debug				= ana.debug;
sM.verbose				= ana.debug;
sM.blend				= true;
sM.bitDepth				= ana.bitDepth;
if exist(ana.gammaTable, 'file')
	clear c
	load(ana.gammaTable);
	if exist('c','var') && isa(c,'calibrateLuminance') && ~isempty(c.finalCLUT) && ~isempty(c.gammaTable)
		sM.gammaTable = c;
	else
		error('CALIBRATION FILE INVALID!!!!!!!!!!!');
	end
	clear c;
	sM.gammaTable.plot;
end
sM.backgroundColour		= ana.backgroundColor;
screenVals				= sM.open; % OPEN THE SCREEN
fprintf('\n--->>> Behavioural Acuity Opened Screen %i : %s\n', sM.win, sM.fullName);
setup(stimuli,sM); %setup our stimulus object

PsychPortAudio('Close');
sM.audio = audioManager(); sM.audio.close();
if IsLinux
	sM.audio.device		= [];
elseif IsWin
	sM.audio.device		= [];
end
sM.audio.setup();
	
%=============================================================
%===========================tobii manager=====================
eT						= tobiiManager();
eT.verbose				= ana.debug;
eT.name					= ana.nameExp;
eT.model				= ana.tracker;
eT.trackingMode			= ana.trackingMode;
eT.eyeUsed				= ana.eyeUsed;
eT.sampleRate			= ana.sampleRate;
eT.calibrationStimulus	= ana.calStim;
eT.manualCalibration	= ana.calManual;
eT.calPositions			= ana.calPos;
eT.valPositions			= ana.valPos;
eT.autoPace				= ana.autoPace;
eT.paceDuration			= ana.paceDuration;
eT.smoothing.nSamples	= ana.nSamples;
eT.smoothing.method		= ana.smoothMethod;
eT.smoothing.window		= ana.w;
eT.smoothing.eyes		= ana.smoothEye;
if ~ana.isDummy; eT.verbose	= true; end
if ~ana.useTracker || ana.isDummy
	eT.isDummy			= true;
end
if length(Screen('Screens')) > 1 && sM.screen - 1 >= 0 && ana.useTracker% ---- second screen for calibration
	s					= screenManager;
	s.verbose			= ana.debug;
	s.screen			= sM.screen - 1;
	s.backgroundColour	= sM.backgroundColour;
	s.distance			= sM.distance;
	[w,h]				= Screen('WindowSize',s.screen);
	s.windowed			= [0 0 round(w/(100/ana.operatorPercent)) round(h/(100/ana.operatorPercent))];
	s.pixelsPerCm		= ana.operatorPPC;
	s.bitDepth			= '';
	s.blend				= sM.blend;
	s.bitDepth			= '8bit';
	s.blend				= true;
end
if exist('s','var')
	initialise(eT,sM,s);
else
	initialise(eT,sM);
end

eT.settings.cal.fixBackSize = round(ana.spotSize * sM.ppd);
eT.settings.cal.fixFrontSize = round(ana.spotLine * sM.ppd);
eT.settings.cal.doRandomPointOrder  = false;

ana.cal=[];
if isempty(ana.calFile) || ~exist(ana.calFile,'file')
	name = regexprep(ana.subject,' ','_');
	ana.calFile = [eT.paths.savedData filesep 'TobiiCal-' name '.mat'];
end
if ana.reloadPreviousCal && exist(ana.calFile,'file')
	load(ana.calFile);
	if isfield(cal,'attempt') && ~isempty(cal.attempt); ana.cal = cal; end
end
cal = trackerSetup(eT, ana.cal); ShowCursor();
if ~isempty(cal) && isfield(cal,'attempt')
	cal.comment=sprintf('Subject:%s | Comments: %s | tobii calibration',ana.subject,ana.comments);
	cal.computer = ana.computer;
	cal.date = ana.date;
	assignin('base','cal',cal); %put our calibration ready to save manually
	if ~exist(ana.calFile,'file')
		save(ana.calFile,'cal');
	end
	ana.outcal = cal;
end

% ---- initial fixation values.
eT.updateFixationValues(ana.XFix,ana.YFix,ana.initTime,ana.fixTime,ana.radius,ana.strict);
resetFixation(eT);

%==========================================================================
%=============================================================TASK Settings
rng('shuffle'); % shuffle the random number generator
if ana.useStaircase == false
	task					= stimulusSequence();
	task.verbose			= ana.debug;
	task.name				= ana.nameExp;
	task.nBlocks			= nBlocks;
	task.nVar(1).name		= 'contrast';
	task.nVar(1).stimulus	= 3;
	task.nVar(1).values		= ana.contrastRange;
	randomiseStimuli(task);
	initialiseTask(task);
else
	task.thisRun			= 0;
	task.minBlocks			= length(ana.contrastRange);
	stopRule				= ana.stopRule;
	task.nRuns				= stopRule;
	usePriors				= ana.usePriors;
	setupStairCase();
end

% ---TASK TYPE
switch ana.myFunc
	case 'Blank Stage 1'
		taskType = 1;
	case 'Blank Stage 2'
		taskType = 2;
	case 'Grating Alone'
		taskType = 3;
	case 'Blank + Grating'
		taskType = 4;
	otherwise
		taskType = 5;
end
ana.taskType = taskType;

%==========================================================================
%=============================================================DRAW EYE
switch lower(ana.drawEye)
	case 'same screen'
		drawEye = 1;
	case 'new screen'
		if exist('s','var') && isa(s,'screenManager')
			if isempty(eT.operatorScreen); eT.operatorScreen = s; eT.secondScreen=true; end
			if ~s.isOpen; s.open; end
			drawEye = 2;
			refRate = 3; %refresh window every N frames
		else
			drawEye = 1;
		end
	otherwise
		drawEye = 0;
end

%=====================================================================
try %our main experimental try catch loop
%=====================================================================
	%===========================prepare===========================
	Priority(MaxPriority(sM.win)); %bump our priority to maximum allowed
	if ~ana.debug; ListenChar(-1); end
	commandwindow;
	startRecording(eT); WaitSecs('YieldSecs',0.1); eT.verbose = false;
	trackerMessage(eT,'!!! Starting Session...');
	breakLoop		= false;
	ana.task		= struct();
	thisRun			= 0;
	breakLoop		= false;
	startAlpha		= stimuli{4}.alphaOut;
	startAlpha2		= stimuli{4}.alpha2Out;
	fadeAmount		= ana.fadeAmount/100;
	fadeAmountTarget = ana.fadeAmountTarget/100;
	fadeFinal		= ana.fadeFinal/100;
	gT				= struct('total',0,'targetTotal',0,'correct',0,...
						'value',0,'run',0,'latest',NaN,...
						'cTotal',zeros(1,task.minBlocks),'cCorrect',zeros(1,task.minBlocks));
	bT				= gT;
	
	%============================================================
	while ~breakLoop
		thisRun							= thisRun + 1;
		%-----setup our values and print some info for the trial
		if thisRun > 1; startRecording( eT ); end
		if ana.useStaircase 
			contrastOut					= staircase.xCurrent;
			ana.task(thisRun).runs		= length(staircase.x);
			ana.task(thisRun).var		= ana.task(thisRun).runs;
			task.totalRuns				= ana.task(thisRun).runs;
		else
			contrastOut					= task.outValues{task.thisRun,1};
			ana.task(thisRun).var		= task.outMap(task.thisRun,1);
			ana.task(thisRun).runs		= [task.totalRuns, task.thisBlock, task.thisRun];
		end
		
		transitionTime = randi( [ana.switchTime(1)*1e3, ana.switchTime(2)*1e3] ) / 1e3;
		targetTime = randi( [ana.targetON(1)*1e3, ana.targetON(2)*1e3] ) / 1e3;
		if taskType == 3
			ana.gProbability			= 100;
		end
		ana.task(thisRun).probability	= rand;
		if ana.task(thisRun).probability <= (ana.gProbability/100) && taskType > 2
			showGrating					= true;
		else
			showGrating					= false;
		end
		
		ana.task(thisRun).thisRun			= thisRun;
		ana.task(thisRun).showGrating		= showGrating;
		ana.task(thisRun).contrast			= contrastOut;
		ana.task(thisRun).transitionTime	= transitionTime;
		ana.task(thisRun).targetTime		= targetTime;
		
		stimuli{2}.xPositionOut			= ana.targetPosition;
		stimuli{3}.contrastOut			= contrastOut;
		stimuli{4}.alphaOut				= startAlpha;
		stimuli{4}.alpha2Out			= startAlpha2;
		stimuli{4}.xPositionOut			= ana.XFix;
		hide(stimuli);
		show(stimuli{4}); % fixation is visi/home/cog5/MatlabFiles/Calibrations/AorusFI27-120Hzble
		update(stimuli);
		
		tStart = 0; tBlank = 0; tGrat = 0; tEnd = 0;
		
		eT.fixInit.X			= [];
		eT.fixInit.Y			= [];
		eT.fixInit.time			= 0.1;
		eT.fixInit.radius		= ana.radius+2;
		eT.updateFixationValues( ana.XFix, ana.YFix, ana.initTime, ana.fixTime, ana.radius, ana.strict );
		resetFixationHistory( eT );
		
		ti='';
		if taskType > 2 && showGrating 
			ti = [ti 'GRATING trial'];	
		elseif taskType > 2 && ~showGrating
			ti = [ti 'BLANK trial'];
		end
		fprintf('\n\n===>>> %s START %i / %i: CONTRAST = %.6f TRANSITION TIME = %.2f TARGET ON = %.2f\n',...
			ti, thisRun, task.totalRuns, contrastOut, transitionTime, targetTime);
		trackerMessage(eT,['TRIALID ' num2str(thisRun)]);
		% ======================================================INITIATE TRIAL
		response			= UNDEFINED; 
		fixated				= '';
		tick				= 1; 
		vbl					= flip(sM); tStart = vbl + sM.screenVals.ifi;
		while ~strcmpi(fixated,'fix') && ~strcmpi(fixated,'breakfix')
			draw(stimuli);
			if drawEye == 1 
				drawEyePosition(eT,true);
			elseif drawEye == 2 && mod(tick,refRate) == 0
				drawGrid(s);
				trackerDrawFixation(eT);
				trackerDrawEyePosition(eT);
			end
			finishDrawing(sM);
		
			vbl = flip(sM);
			if tick == 1; trackerMessage(eT,'INITIATE_FIX',vbl); end
			if drawEye == 2 && mod(tick,refRate)==0; flip(s,[],[],2);end
			tick			= tick + 1;
			
			getSample(eT);
			fixated			= testSearchHoldFixation(eT,'fix','breakfix');
			doBreak			= checkKeys();
			if doBreak; break; end
		end
		if ~strcmpi(fixated,'fix')
			Screen('Flip', sM.win); %flip the buffer
			response		= BREAKINIT;
		else
			trackerMessage(eT,'END_FIX',vbl)
		end
		fprintf('--->>> Time delta Init = %.3f\n',vbl - tStart);
			
		% ======================================================TASKTYPE>0 (Blank)
		if taskType > 0 && response > -1
			if showGrating
				gT.total		= gT.total + 1;
			else
				bT.total		= bT.total + 1;
			end
			tick				= 1;
			triggerTarget		= true;
			triggerFixOFF		= true;
			show(stimuli{1});
			tBlank = vbl + sM.screenVals.ifi; 
			while vbl < (tBlank + transitionTime) && response ~= BREAKBLANK
				
				thisT			= vbl - tBlank;
				
				draw(stimuli); %draw stimulus
				if drawEye == 1 
					drawEyePosition(eT,true);
				elseif drawEye == 2 && mod(tick,refRate) == 0
					drawGrid(s);trackerDrawFixation(eT);trackerDrawEyePosition(eT);
					trackerDrawText(eT,'Blank Period...');
				end
				finishDrawing(sM);
				
				if triggerTarget && thisT > targetTime
					triggerTarget = false; show(stimuli{2}); 
				end
				
				if triggerFixOFF && taskType > 1 && thisT > ana.fixOFF
					stimuli{4}.alphaOut		= stimuli{4}.alphaOut - fadeAmount;
					stimuli{4}.alpha2Out	= stimuli{4}.alpha2Out - fadeAmount;
					if stimuli{4}.alphaOut <= fadeFinal+fadeAmount
						if fadeFinal <= 0
							hide(stimuli{4});
						end
						triggerFixOFF = false;
					end
				end

				vbl = Screen('Flip',sM.win, vbl + screenVals.halfisi); %flip the buffer
				if drawEye == 2 && mod(tick,refRate) == 0; flip(s,[],[],2); end
				tick			= tick + 1;
				getSample(eT);
				isfix			= isFixated(eT);
				if ~isfix
					fixated		= 'breakfix';
					response	= BREAKBLANK;
					statusMessage(eT,'Subject Broke Fixation!');
					trackerMessage(eT,'MSG:BreakFix')
					break;
				end
			end
		end
		if tBlank > 0; fprintf('--->>> Time delta blank = %.3f\n',vbl - tBlank); end
		
		%====================================================TASKTYPE > 2 GRATING/BLANK
		if taskType > 2 && response > -1
			if showGrating
				gT.targetTotal			= gT.targetTotal + 1;
				stimuli{1}.hide(); stimuli{3}.show()
				stimuli{4}.xPositionOut = ana.targetPosition;
				eT.fixInit.X			= ana.XFix;
				eT.fixInit.Y			= ana.YFix;
				eT.updateFixationValues(ana.targetPosition,ana.YFix,ana.initTarget,...
					ana.fixTarget,ana.radius+1,ana.strict);
				fixated					= 'searching';
				fixOFF					= ana.fixOFF2;
				if fixOFF <= 0 %user put 0, means don't show fixation cross over target
					stimuli{4}.alphaOut = 0;
					stimuli{4}.alpha2Out = 0;
					stimuli{4}.hide();
					triggerFixOFF		= false;
				else
					stimuli{4}.alphaOut = startAlpha;
					stimuli{4}.alpha2Out = startAlpha2;
					triggerFixOFF		= true;
				end
				update(stimuli);
			else
				bT.targetTotal			= bT.targetTotal + 1;
				fixOFF					= ana.fixOFF;
				eT.fixation.radius		= ana.radius+2;
				eT.fixation.time		= ana.keepBlank;
				resetFixationTime(eT); fixated = 'fixing';
			end
			
			tGrat						= GetSecs;
			fprintf('--->>> Time delta to switch = %.3f\n',tGrat - vbl);
			while ~strcmpi(fixated,'fix') && ~strcmpi(fixated,'breakfix') && ~strcmpi(fixated,'EXCLUDED!')
				
				if showGrating
					thisT				= vbl - tGrat;
				else
					thisT				= vbl - tBlank;
				end
				
				draw(stimuli); %draw stimulus
				
				if drawEye == 1 
					drawEyePosition(eT,true);
				elseif drawEye == 2 && mod(tick, refRate) == 0
					drawGrid(s);
					trackerDrawFixation(eT);
					trackerDrawEyePosition(eT);
					trackerDrawText(eT, 'Target Period...');
				end
				
				finishDrawing(sM);
				
				if triggerFixOFF && thisT > fixOFF
					stimuli{4}.alphaOut		= stimuli{4}.alphaOut - fadeAmountTarget;
					stimuli{4}.alpha2Out	= stimuli{4}.alpha2Out - fadeAmountTarget;
					if stimuli{4}.alphaOut <= fadeFinal + fadeAmountTarget 
						if fadeFinal <= 0
							hide(stimuli{4});
						end
						triggerFixOFF		= false;
					end
				end
				
				doBreak				= checkKeys();
				if doBreak; break; end
				
				vbl = Screen('Flip', sM.win, vbl + screenVals.halfisi); %flip the buffer
				if drawEye == 2 && mod(tick,refRate) == 0; flip(s,[],[],2);end
				tick				= tick + 1;
				
				getSample(eT);
				if showGrating
					fixated			= testSearchHoldFixation(eT,'fix','breakfix');
				else
					fixated			= testHoldFixation(eT,'fix','breakfix');
				end				
				if strcmpi(fixated,'fix') && showGrating
					response		= YESTARGET;
					if ~ana.useStaircase
						gT.cTotal(ana.task(thisRun).var) = gT.cTotal(ana.task(thisRun).var) + 1;
					end
				elseif strcmpi(fixated,'fix') && ~showGrating
					response		= YESBLANK;
				elseif strcmpi(fixated,'breakfix')
					response		= BREAKTARGET;
					if ~ana.useStaircase
						gT.cTotal(ana.task(thisRun).var) = gT.cTotal(ana.task(thisRun).var) + 1;
					end
				elseif strcmpi(fixated,'EXCLUDED!')
					response		= BREAKEXCL;
					fprintf('---!!! Fix INIT Exclusion triggered!\n');
				end
			end
		end
		tEnd						= flip(sM);
		ana.task(thisRun).tStart	= tStart;
		ana.task(thisRun).tBlank	= tBlank;
		ana.task(thisRun).tGrat		= tGrat;
		ana.task(thisRun).tEnd		= tEnd;
		ana.task(thisRun).RT		= ( tEnd - tGrat ) - ana.fixTarget;
		ana.task(thisRun).xAll		= eT.xAll;
		ana.task(thisRun).yAll		= eT.yAll;
		if tGrat > 0; fprintf('--->>> Time delta grating/blank choice = %.3f RT %.3f | TOTAL = %.3f\n', ...
				tEnd - tGrat, ana.task(thisRun).RT, tEnd - tStart);end
		%====================================================== FINALISE TRIAL
		timeOut						= 1;
		updateResponse();
		
		gT.value					= ( gT.correct / gT.total ) * 100;
		gT.valueTarget				= ( gT.correct / gT.targetTotal ) * 100;
		bT.value					= ( bT.correct / bT.total ) * 100;
		bT.valueTarget				= ( bT.correct / bT.targetTotal ) * 100;
		
		doPlot();
		
		trackerMessage(eT, ['TRIAL_RESULT ' num2str(response)]);
		fprintf('--->>> RESPONSE: %i | Grating Correct: %.2f | Blank Correct: %.2f | TOTAL: %i',...
			response, gT.value, bT.value, thisRun);
		fprintf(' > TRIALS: grating = %i (%.2f) blank = %i (%.2f)\n',...
			gT.total, gT.total / (gT.total + bT.total), bT.total, bT.total / (gT.total + bT.total));
		stopRecording( eT );
		drawBackground( sM );
		vbl							= flip( sM ); %flip the buffer
		t							= vbl;
		while vbl < ( t + timeOut )
			doBreak					= checkKeys();
			if doBreak; break; end
			vbl						= flip( sM );
		end
		
		if ana.useStaircase
			if staircase.stop; breakLoop = true; end
		else
			if task.taskFinished; breakLoop = true; end
		end
		
	end %=====================================================END MAIN WHILE
	
	%=========================================================CLEANUP
	flip(sM);
	Priority(0); ListenChar(0); ShowCursor;
	reset(stimuli); %reset our stimuli
	try if isa(sM,'screenManager'); close(sM); end; end %#ok<*TRYNC>
	try if isa(s,'screenManager'); close(s); end; end
	try if isa(eT,'tobiiManager'); close(eT); end; end
	try if isa(rM,'arduinoManager'); close(rM); end; end
	
	f = figure('Position',[10 10 500 800]);
	t = tiledlayout(f,3,2,'TileSpacing','tight');
	xAll1 = []; xAll2 = []; xAll3 = []; xAll4 = []; xAll5 = [];
	yAll1 = []; yAll2 = []; yAll3 = []; yAll4 = []; yAll5 = [];
	for i = 1:length(ana.task)
		if ana.task(i).response == BREAKINIT
			xAll1 = [xAll1 ana.task(i).xAll];
			yAll1 = [yAll1 ana.task(i).yAll];
		elseif ana.task(i).response == BREAKBLANK
			xAll2 = [xAll2 ana.task(i).xAll];
			yAll2 = [yAll2 ana.task(i).yAll];
		elseif any(ana.task(i).response == [BREAKTARGET BREAKEXCL])
			xAll3 = [xAll3 ana.task(i).xAll];
			yAll3 = [yAll3 ana.task(i).yAll];
		elseif ana.task(i).response == YESBLANK
			xAll4 = [xAll4 ana.task(i).xAll];
			yAll4 = [yAll4 ana.task(i).yAll];
		elseif ana.task(i).response == YESTARGET
			xAll5 = [xAll5 ana.task(i).xAll];
			yAll5 = [yAll5 ana.task(i).yAll];
		end
	end
	tit	= {'BREAKINIT', 'BREAKBLANK', 'BREAKFIX', 'YESBLANK', 'YESTARGET'};
	title(t, 'Eye Positions');
	for i = 1:5
		nexttile( t, i )
		plot(eval(['xAll' num2str(i)]),eval(['yAll' num2str(i)]),'ro');
		xlabel('X Eye Position');
		ylabel('Y Eye Position');
		title( tit{i} );
		xlim([-20 20]); ylim([-20 20]); axis square; grid on; grid minor;
	end
	ana.gT						= gT;
	ana.bT						= bT;
	if ~ana.useStaircase
		ana.staircase = [];
	else
		ana.staircase = staircase;
	end
	if isa(task, 'stimulusSequence')
		ana.response			= task.response;
		ana.responseInfo		= task.responseInfo;
	end
	p = uigetdir(pwd, 'Select Directory to Save Data, CANCEL to NOT SAVE.');
	if ischar(p)
		cd(p);
		save([ana.nameExp '.mat'], 'ana', 'task', 'sM', 'stimuli', 'eT');
		disp(['=====SAVE, saved current data to: ' pwd]);
	else
		eT.saveFile				= ''; %blank save file so it doesn't save
	end	
	assignin('base', 'ana', ana);
catch ME
	getReport(ME)
	try if isa(sM, 'screenManager'); close(sM); end; end %#ok<*TRYNC>
	try if isa(s, 'screenManager'); close(s); end; end
	try if isa(eT, 'tobiiManager'); close(eT); end; end
	try if isa(rM, 'arduinoManager'); close(rM); end; end
	Priority(0); ListenChar(0); ShowCursor;
	disp(['!!!!!!!!=====CRASH, save current data to: ' pwd]);
	save([ana.nameExp 'CRASH.mat'], 'task', 'ana', 'sM', 'stimuli', 'eT', 'ME')
	reset(stimuli);
	clear stimuli sM task taskB taskW md eT s
	rethrow(ME);
end
	
	%=========================================================CHECK KEYS
	function doBreak = checkKeys()
		doBreak					= false;
		[keyIsDown, ~, keyCode] = KbCheck(-1);
		if keyIsDown == 1
			rchar = KbName(keyCode);if iscell(rchar); rchar=rchar{1}; end
			switch lower(rchar)
				case {'q'}
					Screen('DrawText', sM.win, '===>>> EXIT!!!',10,10);
					flip(sM);
					fprintf('===>>> EXIT!\n');
					fixated		= 'breakfix';
					%response	= BREAKINIT;
					breakLoop	= true;
					doBreak		= true;
				case {'p'}
					fprintf('===>>> PAUSED!\n');
					Screen('DrawText', sM.win, '===>>> PAUSED, press key to resume!!!',10,10);
					flip(sM); 
					WaitSecs('YieldSecs',0.15);
					KbWait(-1);
					fixated		= 'breakfix';
					%response	= BREAKINIT;
					doBreak		= true;
				case {'c'}
					WaitSecs('YieldSecs',0.1);
					fprintf('\n\n--->>> Entering calibration mode...\n');
					trackerSetup(eT,eT.calibration);
					fixated		= 'breakfix';
					%response	= BREAKINIT;
					doBreak		= true;
			end
		end
	end

	%=========================================================CHECKREPSONSE
	function updateResponse()
		%===========================================CORRECT BLANK/TARGET
		if response > 0
			if ana.task(thisRun).contrast > 0
				rM.timedTTL();
				sM.audio.beep(1000,0.1,0.1);
			elseif ana.task(thisRun).contrast == 0
				sM.audio.beep(100,0.75,0.75);
			end
			ana.task(thisRun).response = response;
			timeOut				= ana.IFI;
			if response == YESTARGET
				gT.correct		= gT.correct + 1;
				if ana.useStaircase 
					staircase	= PAL_AMPM_updatePM(staircase, true);
				else
					task.updateTask(response,tEnd);
					gT.cCorrect(ana.task(thisRun).var) = gT.cCorrect(ana.task(thisRun).var) + 1;
				end
			elseif response == YESBLANK
				bT.correct		= bT.correct + 1;
			end
			if drawEye == 2
				drawGrid(s);
				trackerDrawFixation(eT);
				trackerDrawEyePositions(eT);
				trackerDrawEyePosition(eT);
				if response == YESTARGET
					trackerDrawText(eT,sprintf('Correct TARGET: %.2f !!!...',gT.value));
				else
					trackerDrawText(eT,sprintf('Correct BLANK: %.2f !!!...',bT.value));
				end
				flip(s,[],[],2);
			end
		%===================================================CORRECT BLANK ONLY
		elseif response == 0 && taskType < 3
			sM.audio.beep(1000,0.1,0.1);
			rM.timedTTL();
			ana.task(thisRun).response	= response;
			timeOut						= ana.IFI;
			bT.correct					= bT.correct + 1;
			if drawEye == 2 
				drawGrid(s);
				trackerDrawFixation(eT);
				trackerDrawEyePositions(eT);
				trackerDrawEyePosition(eT);
				trackerDrawText(eT,'Correct blank!!!...');
				flip(s,[],[],2);
			end
		%=======================================================INCORRECT
		else
			if ana.task(thisRun).contrast == 0 && response == BREAKTARGET
				rM.timedTTL();
				sM.audio.beep(1000,0.1,0.1);
				fprintf('--->>> Reward given as contrast == 0, although trial is incorrect!\n');
			elseif response ~= BREAKINIT
				sM.audio.beep(100,0.75,0.75);
			end
			ana.task(thisRun).response = response;
			timeOut					= ana.punishIFI;
			if (response == BREAKTARGET) 
				if ana.useStaircase
					staircase		= PAL_AMPM_updatePM(staircase, false);
				else
					task.updateTask(response,tEnd);
				end
			end
			if drawEye==2 
				drawGrid(s);
				trackerDrawFixation(eT);
				trackerDrawEyePositions(eT);
				trackerDrawEyePosition(eT);
				if response == BREAKINIT
					trackerDrawText(eT,'Break in INIT!!!');
					timeOut = 0.75;
				elseif response == BREAKBLANK
					trackerDrawText(eT,'Break in BLANK!!!');
				elseif response == BREAKTARGET
					trackerDrawText(eT,'Break in TARGET!!!');
				elseif response == BREAKEXCL
					trackerDrawText(eT,'Break in TARGET EXCLUSION!!!');
				end
				flip(s,[],[],2);
			end
		end
	end
	
	%=========================================================PLOT
	function doPlot()
		
		cla(ana.plotAxis1);
		hold(ana.plotAxis1,'on')
		if gT.total > 0
			gT.run(gT.total) = gT.value; 
			try plot(ana.plotAxis1, 1:length(gT.run), gT.run,'LineWidth',1,'Color',[0 0.4470 0.7410]); end
			if gT.targetTotal > 0
				gT.runTotal(gT.targetTotal) = gT.valueTarget; 
				try plot(ana.plotAxis1, 1:length(gT.runTotal), gT.runTotal,'--o','LineWidth',1,'Color',[0 0.4470 0.7410]); end
				if gT.targetTotal > 6
					gT.latest = mean(gT.runTotal(end-6:end));
				end
			end

		end
		if bT.total > 0
			bT.run(bT.total) = bT.value; 
			try plot(ana.plotAxis1, 1:length(bT.run), bT.run,'LineWidth',1,'Color',[0.8500 0.3250 0.0980]); end
			if bT.targetTotal > 0
				bT.runTotal(bT.targetTotal) = bT.valueTarget;
				try plot(ana.plotAxis1, 1:length(bT.runTotal), bT.runTotal,'--p','LineWidth',1,'Color',[0.8500 0.3250 0.0980]); end
				if bT.targetTotal > 6
					bT.latest = mean(bT.runTotal(end-6:end));
				end
			end
		end
		hold(ana.plotAxis1,'off');
		box(ana.plotAxis1, 'on'); grid(ana.plotAxis1, 'on');
		xlabel(ana.plotAxis1,'Runs');
		ylim(ana.plotAxis1,[-5 105]);
		ylabel(ana.plotAxis1,'Correct %');
		title(ana.plotAxis1,sprintf('Last 6 Grat: %.2f | Last 6 Blank: %.2f',gT.latest,bT.latest));
		list = {};
		if gT.total > 0;		list{end+1} = 'Grating';	end
		if gT.targetTotal > 0;	list{end+1} = 'Grating *';	end
		if bT.total > 0;		list{end+1} = 'Blank';		end
		if bT.targetTotal > 0;	list{end+1} = 'Blank *';	end
		legend(ana.plotAxis1,list,'Location','southwest');
		
		if ana.useStaircase == true
			cla(ana.plotAxis2); hold(ana.plotAxis2,'on');
			if ~isempty(staircase.threshold)
				rB = linspace(min(staircase.stimRange),max(staircase.stimRange),200);
				outB = ana.PF([staircase.threshold(end) ...
					staircase.slope(end) staircase.guess(end) ...
					staircase.lapse(end)], rB);
				plot(ana.plotAxis2,rB,outB,'r-','LineWidth',2);
				r = staircase.response;
				t = staircase.x(1:length(r));
				yes = r == 1;
				no = r == 0; 
				plot(ana.plotAxis2,t(yes), ones(1,sum(yes)),'ko','MarkerFaceColor','r','MarkerSize',10);
				plot(ana.plotAxis2,t(no), zeros(1,sum(no))+ana.gamma,'ro','MarkerFaceColor','w','MarkerSize',10);
				box(ana.plotAxis2, 'on'); grid(ana.plotAxis2, 'on');
				ylim(ana.plotAxis2, [ana.gamma-0.05 1+0.05]);
				xlim(ana.plotAxis2, [0 max(rB)]);
				xlabel(ana.plotAxis2, 'Contrast');
				ylabel(ana.plotAxis2, 'Responses');
				hold(ana.plotAxis2, 'off');		
			end
		else
			cla(ana.plotAxis2);
			x = task.nVar(1).values;
			y = (gT.cCorrect ./ gT.cTotal) * 100;
			y(isnan(y)) = 0;
			hold(ana.plotAxis2, 'on');	
			plot(ana.plotAxis2, x,y, 'k-o');
			ylim(ana.plotAxis2, [-5 105]);
			title(ana.plotAxis2, 'Contrast Responses')
			xlabel(ana.plotAxis2, 'Contrast'); 
			ylabel(ana.plotAxis2, '% Correct');
			hold(ana.plotAxis2, 'off');	
		end
		drawnow;
	end

	function setupStairCase()
		priorAlpha = linspace(min(ana.contrastRange), max(ana.contrastRange),ana.alphaGrain);
		priorBeta = linspace(0, ana.betaMax, ana.betaGrain); %our slope
		priorGammaRange = ana.gamma;  %fixed value (using vector here would make it a free parameter)
		priorLambdaRange = ana.lambda; %ditto
		
		staircase = PAL_AMPM_setupPM('stimRange',ana.contrastRange,'PF',ana.PF,...
			'priorAlphaRange', priorAlpha, 'priorBetaRange', priorBeta,...
			'priorGammaRange',priorGammaRange, 'priorLambdaRange',priorLambdaRange,...
			'numTrials', stopRule,'marginalize',ana.marginalize);
		
		if usePriors
			prior = PAL_pdfNormal(staircase.priorAlphas,ana.alphaPrior,ana.alphaSD).*PAL_pdfNormal(staircase.priorBetas,ana.betaPrior,ana.betaSD);
			figure;
			imagesc(staircase.priorBetaRange,staircase.priorAlphaRange,prior);axis square
			ylabel('Threshold');xlabel('Slope');title('Initial Bayesian Priors')
			staircase = PAL_AMPM_setupPM(staircase,'prior',prior);
		end
	end

	function saveMetaData()
		ana.values.nBlocks			= nBlocks;
		ana.values.nBlocksOverall	= nBlocksOverall;
		ana.values.NOSEE			= YESBLANK;
	    ana.values.YESSEE			= YESTARGET;
		ana.values.UNSURE			= UNSURE;
		ana.values.BREAKFIX			= BREAKTARGET;
		ana.values.BREAKFIXEXCL		= BREAKEXCL;
		ana.values.BREAKINIT		= BREAKINIT;
		ana.values.BREAKBLANK		= BREAKBLANK;
		ana.values.UNDEFINED		= UNDEFINED;
	end

end


