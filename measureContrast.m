%% This file runs a grating and then measures the values from the screen as the grating is moved

gratingType = 'colour';
contrast = 0.25;
bitDepth = 'EnableBits++Color++Output'; %for display++
%bitDepth = 'Native10bit'; %for aorus
comment = 'Aorus in 10bit mode';
needGamma = true;
background = [0.25 0.25 0.25];
colour = [1 1 1];
colour2 = [0 0 0];
resolution = 2^10;
testSpot = false;
flipTest = true;

s = screenManager;
d = discStimulus;
if strcmpi(gratingType,'colour')
	g = colourGratingStimulus;
else
	g = gratingStimulus;
end
if needGamma
	%load('~/MatlabFiles/Calibrations/TobiiTX300_SET2_MonitorCalibration.mat');
	load('~/Code/Training/AorusFI27-120Hz-NEWcalibration.mat');
	%load('~/MatlabFiles/Calibration/Display++Color++Mode-Ubuntu-RadeonPsychlab.mat')
	%c.choice = 2;
	%c.plot;
	s.gammaTable = c;
else
	c = calibrateLuminance;
end

% setup screen
s.bitDepth = bitDepth;
s.backgroundColour = background;
sv = open(s);

%g.type='square';
g.mask = false;
g.size = 30;
g.tf=0;
g.sf = 0.2;
g.contrast = contrast;
if isprop(g,'colour2')
	g.colour = colour;
	g.cololur2 = colour2;
	g.correctBaseColour = true;
else
	g.colour = s.backgroundColour;
end
d.colour = [0.2 0.2 0.2];
d.size = 30;

% setup the stimuli
c.screenVals = sv;
setup(g, s);
setup(d, s);
draw(g);
flip(s);

% open spectrocal
c.openSpectroCAL();
c.spectroCalLaser(true);
input('Align Laser then press enter to start...')
c.spectroCalLaser(false);

Priority(MaxPriority(s.win));
WaitSecs(0.5);

if flipTest
	g.sf = 0.1; g.update;
	ctest = [0.0001 0.0005 0.001];
	clear Y YY YYY A B;
	h=figure;
	tl = tiledlayout(h,2,length(ctest));
	mn = inf;
	mx = -inf;
	for loop = 1:length(ctest)
		
		g.driftPhase = 0;
		g.contrastOut = ctest(loop);
		
		for i = 1:10
			g.draw();
			s.flip();
			WaitSecs(0.1);
			[~,~,Y(i)] = c.getSpectroCALValues();
			g.driftPhase = g.driftPhase + 180;
		end
	
		A = Y(1:2:9);
		B = Y(2:2:10);
		ct(loop).A = A;
		ct(loop).B = B;
		
		if mn > min([min(A) min(B)])
			mn = min([min(A) min(B)]);
		end
		if mx < max([max(A) max(B)])
			mx = max([max(A) max(B)]);
		end
	
		nexttile(loop);
		plot(A,'ko');
		hold on;
		plot(B,'ro');
		box on;grid on
		title(['Contrast: ' num2str(ctest(loop))]);
		nexttile(loop+length(ctest));
		boxplot([A,B],[ones(1,5),ones(1,5)*2],'Notch','on','Labels',{'phase0','phase180'});
		box on;grid on
		title(['Contrast: ' num2str(ctest(loop))]);
		drawnow;
	end
	
	for loop = 1:length(ctest)
		ax = nexttile(loop);
		ax.YLim = [mn mx];
		ax = nexttile(loop+length(ctest));
		ax.YLim = [mn mx];
	end
	
	tl.YLabel.String = 'Luminance (cd/m^2)';
	tl.Title.String = 'Contrasts:';
	
	for i = 1:length(ct)-1	
		[~,p] = ttest2(ct(i).A,ct(i+1).A);
		[~,p2] = ttest2(ct(i).B,ct(i+1).B);
		fprintf('\nLOW = %.2f : %.2f  p = %.4f | HIGH = %.2f : %.2f p = %.4f\n',...
			median(ct(i).A),median(ct(i+1).A),p,...
			median(ct(i).B),median(ct(i+1).B),p2);
	end
end

phs = [0:22.5:360];
YY=[];
g.contrastOut = 0.005;
g.driftPhase = phs(1);
g.draw();
s.flip();
WaitSecs(0.5);

for loop = 1:length(phs)
	g.driftPhase = phs(loop);
	g.draw();
	s.flip();
	WaitSecs(0.2);
	[~,~,YY(loop)] = c.getSpectroCALValues();
	fprintf('Phase is: %.2f, Luminance is %.4f\n',phs(loop),YY(loop));
end

h=figure('Name',comment);
tl = tiledlayout(h,'flow');
nexttile
plot(phs,YY,'r-o');box on;grid on
title(['Contrast: ' num2str(g.contrastOut)]);
xlabel('Phase (deg)')
ylabel('Output Luminance (cd/m^2)');
drawnow;

s.flip();
range = 0:1/resolution:1;
steps = floor(length(range)/2):floor(length(range)/2)+25;
for loop = steps
	d.colourOut = [range(loop) range(loop) range(loop) 1];
	d.update();
	d.draw();
	s.flip();
	WaitSecs(0.2);
	[~,~,YYY(loop)] = c.getSpectroCALValues();
	fprintf('Loop %i - In/out Luminance %.4f = %.4f\n',loop,range(loop),YYY(loop));
end
YYY=YYY(steps);
nexttile
plot(range(steps),YYY,'r-o');
title([s.bitDepth ' Luminance: ' num2str(resolution)]);
xlabel('Grayscale Step 0-1')
ylabel('Output Luminance (cd/m^2)');
box on;grid on
drawnow;

ListenChar(0);ShowCursor;Priority(0);
g.reset;
d.reset;
s.close;
c.closeSpectroCAL
c.close;
