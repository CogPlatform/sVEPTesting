function [trl, events, triggers] = loadCOGEEG(cfg)
% Uses 8 analog recorded channels to encode strobed words to generate trial
% markers

% read the header information and the events from the data
if isfield(cfg,'header')
	hdr = cfg.header;
else
	hdr   = ft_read_header(cfg.dataset);
end
% read the events from the data
if ~isfield(cfg,'denoise')
	cfg.denoise = true;
end
chanindx      = cfg.chanindx;
detectflank   = 'up';
threshold     = cfg.threshold; % or, e.g., 1/2 times the median for down flanks
event = ft_read_event(cfg.dataset,'header',hdr,...
		'detectflank',detectflank,'chanindx',chanindx,...
		'threshold',threshold,'denoise',cfg.denoise);
labels = hdr.label(chanindx);
nChannels = length(labels);
events = [];
time = linspace(0, (1/hdr.Fs)*hdr.nSamples, hdr.nSamples);

% any trigger < minNextTrigger ms after previous is considered artifact and removed
minNextTrigger = cfg.minTrigger * 1e-3;
% number of samples to allow jitter to assign to same strobed word
nSamples = cfg.jitter; 
% pre stimulus time to select
preTime = cfg.preTime;

%parse our events, removing any events < minNextTriggerTime
for i = 1:nChannels
	lidx = cellfun(@(x)strcmpi(x,labels{i}),{event.type});
	events(i).label		= labels{i};
	events(i).idx		= find(lidx==1);
	events(i).evnt		= event(events(i).idx);
	events(i).samples	= [events(i).evnt.sample];
	events(i).times		= time(events(i).samples);
	events(i).rmIdx		= [];
	times = events(i).times;
	if ~any(diff(times) <= minNextTrigger);continue;end
	j = 1;
	while j < length(times)
		% this checks if the next time < min and also next+1
		% as some TTL signals show ringing.
		if (times(j+1)-times(j)) <= minNextTrigger
			if (j<length(times)-1) && (times(j+2)-times(j)) <= minNextTrigger
				events(i).rmIdx = [events(i).rmIdx j+1 j+2];
				j = j + 3;
			else
				events(i).rmIdx = [events(i).rmIdx j+1];
				j = j + 2;
			end
		else
			j = j + 1;
		end
	end
	if ~isempty(events(i).rmIdx)
		events(i).idx(events(i).rmIdx) = [];
		events(i).evnt(events(i).rmIdx) = [];
		events(i).samples(events(i).rmIdx) = [];
		events(i).times(events(i).rmIdx) = [];
	end
end

% make strobed words from individual events
triggers = [];
times = [];
bidx = 1;

for i = 1:nChannels
	for j = 1:length(events(i).idx)
		fixit = zeros(1,length(events));
		if ~any(times == events(i).times(j)) % check our time list of previously measured events
			triggers(bidx).time = events(i).times(j);
			triggers(bidx).sample = events(i).samples(j);
			triggers(bidx).bword = '00000000';
			triggers(bidx).bword(i) = '1';
			fixit(i) = j;
			% for each event, now we check all other channels for events
			% within 4ms
			for k =  1 : length(events) %check all other channels
				if k == i; continue; end
				[idx,val,delta] = findNearest(events(k).samples,triggers(bidx).sample);
				if delta <= nSamples
					triggers(bidx).bword(k) = '1';
					triggers(bidx).time = min([val,triggers(bidx).time]);
					triggers(bidx).sample = min([events(i).samples(j),events(k).samples(idx)]);
					fixit(k) = idx;
				end
			end
			% make sure all other channels use the same time and sample so
			% we don't double count events
			for l = 1:length(fixit)
				if fixit(l) > 0
					events(l).times(fixit(l)) = triggers(bidx).time;
					events(l).samples(fixit(l)) = triggers(bidx).sample;
				end
			end
			times(end+1) = triggers(bidx).time;
			triggers(bidx).bword = fliplr(triggers(bidx).bword);
			triggers(bidx).value = bin2dec(triggers(bidx).bword);
			triggers(bidx).map = [i j];
			triggers(bidx).fixit = fixit;
			bidx = bidx + 1;
		end
	end
end
%TODO 0 triggers
% resort by times
trl = [];
if isempty(triggers)
	warning('No triggers found');
	trl = [1 hdr.nSamples];
	return
end
[~,sidx] = sort([triggers.time]);
triggers = triggers(sidx);
% make sure any events with identical times are removed
purgeIdx = find(diff([triggers.time]) == 0)+1; 
if ~isempty(purgeIdx)
	warning('Duplicate events! removing...')
	triggers(purgeIdx) = [];
end

% now we need to make the trl structure fieldtrip needs:
% first find all trials where a number is followed by 255
bSamples = round(preTime / (1/hdr.Fs));
nTriggers = length(triggers);
trlN = 0;
trl = [];
for i = 1:(nTriggers - 1)
	if isempty(cfg.correctID) || cfg.correctID==0
		if triggers(i).value ~= 255 && triggers(i+1).value == 255
			trlN = trlN + 1;
			trl(trlN,1) = triggers(i).sample-bSamples;
			trl(trlN,2) = triggers(i+1).sample+bSamples;
			trl(trlN,3) = -bSamples;
			trl(trlN,4) = triggers(i).value;
		end
	else
		if triggers(i).value ~= 255 && triggers(i+1).value == 255 && triggers(i+2).value == cfg.correctID
			trlN = trlN + 1;
			trl(trlN,1) = triggers(i).sample-bSamples;
			trl(trlN,2) = triggers(i+1).sample+bSamples;
			trl(trlN,3) = -bSamples;
			trl(trlN,4) = triggers(i).value;
		end
	end
end

% now we remove duplicate numbers (early bug) or any > 10
%trlN = 0;
%trl = [];
%for i = 1:size(trl0,1)-1
%	if (trl0(i,4) ~= trl0(i+1,4)) && trl0(i,4) <= 249
%		trlN = trlN + 1;
%		trl(trlN,:) = trl0(i,:);
%	end
%end

function [idx,val,delta]=findNearest(in,value)
	%find nearest value in a vector, if more than 1 index return the first	
	[~,idx] = min(abs(in - value));
	val = in(idx);
	delta = abs(value - val);
end


end
