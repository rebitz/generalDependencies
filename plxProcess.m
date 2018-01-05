function [trials,t] = plxProcess(plxfile,trials,startEvent,stopEvent)
%%
% trials = plxProcess(plxfile,trials,startEvent,stopEvent)
%
% plxProcess
%   - append spikes (from file "plxfile") to data struct "trials"
%
% TO DO: add TRACKER TIME OFFSET marker - diff btw plex and ctx, plus scale
%   - to aid in aligning spikes to CTX time events (like saccades)
%
% -- add some constraints around the fitting -> i.e. first/second half of
% plx file for the MGS data

pad = 0.5; % in s, time on either side of start/stop

noStops = 0;

if nargin < 3
    % set some parameters
    startEvent = 'EVT01'; % start FHC
    stopEvent = 'EVT05'; % stop FHC
elseif nargin < 4
    stopEvent = 'EVTdummy'; % stop FHC
    disp('no stop given, using start-to-start time!')
    noStops = 1;
end

if isfield(trials,'trialstart')
    startfield = 'trialstart';
elseif isfield(trials,'trialStart')
    startfield = 'trialStart';
end

% eventNames = {'juice','trialStart','fixOn','fixAcq','targOn',... % 5
%     'targChange','targAcq','trialStop'}; % these will appear in order

if strcmp(stopEvent,'EVT04');
    stopfield = 'fixAcq';
elseif strcmp(stopEvent,'EVT07');
    stopfield = 'targAcq';
elseif ~strcmp(stopEvent,'EVTdummy');
    if isfield(trials,'trialstop')
        startfield = 'trialstop';
    elseif isfield(trials,'trialStop')
        stopfield = 'trialStop';
    else
        disp('no stop field!')
        return
    end
end

%%
% other event information - for reference
% 1. trial start
% 2. juice
% 3. stimulation (grass)
% 4. stimulation (ctx)
% 5. trial stop
% 6. fixation
% 7. target
% 8. go cue
% 9. correct trial marker

% echo the plexon file name
plxfile

% from plexon: channel n is at index n+1
[tscounts,junk,evtcounts] = plx_info(plxfile,1);
chIdx = find(sum(tscounts)>0); chIdx = chIdx - 1;
evtIdx = find(evtcounts > 0);

% get all the events and put them into a struct
t = ([]);
for i = 1:length(evtIdx)
    [n,ts] = plx_event_ts(plxfile,evtIdx(i));
    if n > 0
        evtName = strcat('EVT0',num2str(evtIdx(i)));
        t.(evtName) = ts;
    end
end

if noStops == 1 % if no stop evt info provided, use the next start
    tmp = [t.(startEvent)];
    t.EVTdummy = [tmp(2:end)];
end

% fix any offset in size
if length([t.(startEvent)]) < length([t.(stopEvent)])
    disp('fixing size offset part 1')
    stops = NaN(1,length([t.(startEvent)]));
    for i = 1:length([t.(startEvent)])-1
        idx = and([t.(stopEvent)]>[t.(startEvent)(i)],[t.(stopEvent)]<[t.(startEvent)(i+1)]);
        if sum(idx)>0
            stops(i) = find(idx,1,'first');
        end
    end
    stops = stops(~isnan(stops));
    [t.(stopEvent)] = [t.(stopEvent)(stops)];
end

if length([t.(stopEvent)]) < length([t.(startEvent)])
    disp('fixing size offset part 2')
    disp('(necessary for a no-stop marker parse)')
    % make a mask for the start events that don't have a stop
    mask = zeros(length([t.(startEvent)]),1);
    for i = 1:length([t.(stopEvent)])
        idx = find([t.(startEvent)]<[t.(stopEvent)(i)],1,'last');
        mask(idx) = 1;
    end
    [t.(startEvent)] = [t.(startEvent)(logical(mask))];
end
% 
% keyboard();
% now figure out the mapping to particular trials
theseTs = [t.(stopEvent)] - [t.(startEvent)]; %theseTs = theseTs*1000; % for ctx
if ~noStops
    goalTs = [trials.(stopfield)]' - [trials.(startfield)]';
else
%     tmp = [trials.(startfield)]';
%     goalTs = [tmp(2:end); NaN] - [trials.(startfield)]';
    theseTs = diff([t.(startEvent)]); %theseTs = theseTs*1000; % for ctx
    goalTs = diff([trials.(startfield)]');
end

% if ~noStops
%     % remove outliers in theseTs just in case 
%     [theseTs(theseTs > nanmean(theseTs)+3*std(theseTs))] = deal(NaN)
% end

if length(goalTs) > length(theseTs)
    fprintf('\n ERROR: plexon file is not long enough \n')
    return
end

i = 1; error = []; % find the alignment that minimizes errror
while i+length(goalTs)-1<=length(theseTs); % while goals still fit in realities
    pred = theseTs(i:i+length(goalTs)-1);
    temp = goalTs - pred;
    error(i) = sum(temp(~isnan(temp)).^2);
    i = i+1;
end

% assign the start and end trial idx's to these points
[junk,startIdx] = min(error);
endIdx = startIdx+length(goalTs)-1;

% evaluate and report on the quality of the alignment
%[h,p] = ttest2(min(error),error);
[p,h] = signtest(error,min(error));
if p < 0.0001
    fprintf('\n Plexon alignment OK, p < %d \n',p)
else
    fprintf('\n CAUTION!!! \n')
    fprintf('\n CAUTION!!! \n')
    fprintf('\n CAUTION!!! \n')
    fprintf('\n P value at %d \n',p)
    fprintf('\n Poor quality plexon alignment - check your data! \n')
    fprintf('\n CAUTION!!! \n')
    fprintf('\n CAUTION!!! \n')
    fprintf('\n CAUTION!!! \n')
    figure(99); plot(gsmooth(error,5)); hold on
    temp = ylim;
    line([startIdx startIdx],[temp(1) temp(2)]);
end

%%
% make a splitter of the data file in case it's broken up somewhere?

% to convert from cortex to plexon timestamps

% initialize event field names
evts = fieldnames(t);
for i = 1:length(evts)
    [trials.(evts{i})] = deal(NaN);
end

% start filling in the trial struct
for i = 1:length(trials);
    
    start = [t.(startEvent)(startIdx+i-1)];
    stop = [t.(stopEvent)(startIdx+i-1)];
    
    trials(i).(startEvent) = start;
    trials(i).(stopEvent) = stop;
    
    % then add all the non-start/stop events
    evts = evts(~strcmp(evts,startEvent));
    evts = evts(~strcmp(evts,stopEvent));
    
    for j = 1:length(evts);
        temp = find(and(t.(evts{j})>start,t.(evts{j})<stop));
        if ~isempty(temp) && length(temp) == 1 % require that it find least 1 evt per trial
            trials(i).(evts{j}) = t.(evts{j})(temp);
        end
    end
end

fprintf('\n Events appended to data struct, getting spike data... \n')

% add spikes to the trial struct
% supply names for the units we'll pull out
unitDict = {'i','a','b','c','d','e','f','g','h'};
zpad = '00'; units = 0;

for i = 1:length(chIdx)

    % knock zero padding down if we go over 10 units
    if chIdx(i) > 9; zpad = '0'; end

    % find all available units on this channel
    unitIdx = tscounts(:,chIdx+1);

    for j = 1:sum(unitIdx>1)

        [n,ts] = plx_ts(plxfile,chIdx(i),j-1);        
        unitName = strcat('sig',zpad,num2str(chIdx(i)),unitDict{j});

        units = units+1;
        
        for k = 1:length(trials)
            start = trials(k).(startEvent)-pad;
            stop = trials(k).(stopEvent)+pad;
            idx = and(ts>start,ts<stop);
            if sum(idx)>0
                trials(k).(unitName) = ts(idx);
            else
                trials(k).(unitName) = NaN;
            end
        end
    end

end

fprintf('\n Spikes appended to data struct, %d signals found. \n', units)

%% CAUTION: method below is NOT GOING TO WORK with multiple channels

if length(chIdx) == 1
    fields = fieldnames(trials);
    idx = and(~cellfun(@isempty,strfind(fields,num2str(chIdx))),...
        ~cellfun(@isempty,strfind(fields,'sig')));
    siglist = fields(idx);
    
    zpad = '00';
    if chIdx > 9; zpad = '0'; end
    newfield = strcat('sig',zpad,num2str(chIdx),'all');
    
    [trials.(newfield)] = deal(NaN);

    for i = 1:length(trials);
        for j = 1:length(siglist)

            sigOI = siglist{j};
            trials(i).(newfield) = [trials(i).(newfield); trials(i).(sigOI)];

        end
    end
end
