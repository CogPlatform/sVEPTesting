%==========================================================================
%===========================================================response values
YESBLANK = 1; YESTARGET = 2; UNSURE = 4; BREAKINIT = -100; BREAKBLANK = -10;
BREAKTARGET = -1; BREAKEXCL = -5; UNDEFINED = 0;

% parametric standard error?
parametric = NaN;

%do we include the break exclude trials?
useExclusion = false;

% do we treat max contrast as only useful for lapse rate, 'nAPLE' (no) or
% 'jAPLE' (yes)
lapseFits = 'nAPLE';

% ---- which psychometric function to use?
PF = @PAL_Gumbel;

% ---- load data
file = uigetfile;
if file == 0; return; end
load(file);
fprintf('\n\DATA: %s\n',file);

% ---- extract values from data
contrasts = task.nVar(1).values;
trials = ana.task([ana.task.showGrating]==true);
trialscorrect = trials([trials.response]==YESTARGET);
trialswrong = trials([trials.response]==BREAKTARGET);
trialswrongall = trials([trials.response]==BREAKTARGET | [trials.response]==BREAKEXCL);

ti=sprintf('Trials: %i - Corr: %i - BREAK: %i - BREAKALL: %i - EXCL: %i',...
	length(trials),length(trialscorrect),length(trialswrong),length(trialswrongall),...
	useExclusion);
disp(ti);

if useExclusion
	trialswrong = trialswrongall;
end

contrastTotal = [];
contrastCorrect = [];

for i = 1 : length(contrasts)
	tr = trialscorrect([trialscorrect.contrast] == contrasts(i));
	contrastCorrect(i) = length(tr);
	tr = trialswrong([trialswrong.contrast] == contrasts(i));
	contrastWrong(i) = length(tr);
	contrastTotal(i) = contrastWrong(i) + contrastCorrect(i);
end

contrasts(contrasts==0) = 1e-6;

total = contrastTotal;
correct = contrastCorrect;

% ================= Here are our model parameters =========================
% ---- threshold
search.alpha = logspace(min(contrasts), max(contrasts), 100);

% ---- slope
search.beta = logspace(0, 2, 100);

% ---- guess rate
search.gamma = 0;
%searchGrid.gamma = linspace(0,0.5,10);

% ---- lapse bias
%search.lambda = 0;
search.lambda = linspace(0,0.2,10);

% ---- which parameters to search
freeParameters = [1 1 0 1];

% ============================ Here we run the model ======================
[params,b,c,d] = PAL_PFML_Fit(contrasts,correct,total,search,freeParameters,PF,...
	'lapseLimits',[0 0.2],'guessLimits',[0 0.5],'lapseFits',lapseFits);

if c ~= 1
	warndlg('DID NOT CONVERGE!!!!!');
end
if params(2) == Inf
	warning('Had to change the INF slope parameter!!')
	params(2) = max(search.beta);
end
xrange = linspace(0, max(contrasts), 1000);
fit = PF(params,xrange);

% ========================= And plot our result ===========================
fname = functions(PF);
fname = fname.function;
r = groot; ss = r.ScreenSize;
f = figure('Position',[ss(3)/4+randi(50) ss(4)/2 ss(3)/2 ss(4)/2]);
tiledlayout(f,'flow');
nexttile; hold on
plot(contrasts,(correct./total),'k.','Color',[0.3 0.3 0.3],'MarkerSize',30);
plot(xrange,fit,'Color',[0.8 0.5 0],'LineWidth', 2);
title([fname ' PARAMS: ' num2str(params,'%.4f ') ' LL: ' num2str(b,'%.3f')],'Interpreter','none');
subtitle([file ' ' ti],'Interpreter','none');
xlim([-(max(contrasts)/100) max(contrasts)+(max(contrasts)/100)]);
ylim([-0.05 1.05]);
xlabel('Stimulus Contrast');
ylabel('Proportion correct [0 - 1]');
grid on; grid minor; box on; hold off
drawnow;

% =======================Get errors and goodness of fit? =================
if islogical(parametric)
	warning off
	t = text(0.01, 0.1, 'Calculating GOF, please wait...','FontSize',16);drawnow
	B=400;
	if parametric == 1
		[SD, paramsSim, LLSim, converged] = PAL_PFML_BootstrapParametric(...
			contrasts, total, params, freeParameters, B, PF, ...
			'searchGrid', search);
	else
		[SD, paramsSim, LLSim, converged] = PAL_PFML_BootstrapNonParametric(...
			contrasts, correct, total, [], freeParameters, B, PF,...
			'searchGrid',search);
	end

	message = sprintf('Threshold SE: %6.4f',SD(1));
	message = [message sprintf(' | slope SE: %6.4f',SD(2))];
	message = [message sprintf(' | lapse SE: %6.4f',SD(4))];
	
	%Number of simulations to perform to determine Goodness-of-Fit
	B=500;
	disp('Determining Goodness-of-fit.....');

	[Dev, pDev] = PAL_PFML_GoodnessOfFit(contrasts, correct, total, ...
		params, freeParameters, B, PF, 'searchGrid', search);

	%Put summary of results on screen
	message = [message sprintf(' | deviance: %6.4f',Dev)];
	message = [message sprintf(' | p: %6.4f',pDev)];
	
	delete(t);
	text(0.01, 0.1, message,'FontSize',14);
	warning on
	
	%Create simple plot
% 	ProportionCorrectObserved = correct./total; 
% 	StimLevelsFineGrain = linspace(min(contrasts), max(contrasts), 500);
% 	ProportionCorrectModel = PF(params,StimLevelsFineGrain);
% 	
% 	nexttile;
% 	hold on;
% 	plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0.8 0.5 0],'linewidth',2);
% 	plot(contrasts,ProportionCorrectObserved,'k.','markersize',30);
% 	axis([min(StimLevels) max(StimLevels) -0.05 1.05]);
% 	xlabel('Stimulus Contrast');
% 	ylabel('Proportion correct [0 - 1]');
% 	title(message)
% 	grid on; grid minor; box on; hold off
end
