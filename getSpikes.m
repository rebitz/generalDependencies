function [meanResponse,stdResponse,meanTrace,steTrace,trialData] = getSpikes(trials,evtOI,sigOI,pre,post)
%
% fast FR analysis and PSTH generation (sparse arrays and logical indexing)
%   processes an entire 2000 trial experiment in < 250 ms
%
% takes inputs
%    1. data struct - with all time units already in ms
%    2. event timestamp pointer (str, field in data struct)
%    3. signal pointer (str, field in data struct)
%    4. prespace (default: 0ms)
%    5. postspace (default: 500ms)
% returns
%    1. average response of that signal in that epoch (spikes/s)
%    2. standard deviation of the response across the epoch
%    3. mean FR trace (PSTH)
%    4. ste of the PSTH
%    5. trial-by-trial vector of FR mean across specified epoch
%
% EXAMPLE: response = getSpikes(data,'evtOI','sigOI',300,0)
%    looks for 'evtOI' and 'sigOI' in data struct 'data'
%    returns the avg response over the 300ms-0ms before the event
%
%
% 12/08/13: written as function, RBE

if nargin < 5
    disp('WARNING:')
    disp('Missing epoch arguments, using the 0-500ms post event.')
    pre = 0;
    post = 500;
end


% if ~isnan(trials.(sigOI)) && range(trials(1).(sigOI)) < 100
%     fprintf('\n CAUTIONS!!! \n')
%     fprintf('\n Suspect time stamps are in seconds, correcting for this. \n')
    scaleby = 1000;
% else
%     scaleby = 1;
% end

% initialize spike pointers for sparse array
x = []; y = [];
xpos = [-pre:1:post];

% go through all selected trials
for i = 1:length(trials);

    % pull out the event timestamp
    tstamp = [trials(i).(evtOI)]*scaleby;

    if ~isempty(tstamp)
        % rescale to ms
        fu = [trials(i).(sigOI)]*scaleby;
        % center them in 1 ms bins
        fu = round(fu);
        % find relevant stamps
        fu = fu(and(fu > (tstamp-pre),fu < (tstamp+post)));
        % center it around the event
        fu = fu - round(tstamp);
        % add the prespace to get rid of -idx
        fu = fu+pre+1;

        y = [y; ones(length(fu),1)*i]; % trial number
        x = [x; fu]; % time stamp of event
    end
    
end

A = sparse(y, x, [ones(length(y),1)], length(trials), length(xpos));
A = A*scaleby; % convert from spikes/ms to spikes/s

meanTrace = full(mean(A,1));

trialData = full(mean(A,2));
meanResponse = mean(full(mean(A,1)));


if length(trials) > 1
%     stdResponse = mean(full(std(A,1))./sqrt(length(trials)-1));
    stdResponse = std(full(mean(A,1))); % changed 10/18/14, was wrong
    steTrace = full(std(A,1)./sqrt(length(trials)-1));
else
%     disp('WARNING:')
%     disp('Too few trials (<1) for a meaningful variance estimate.')
%     disp('Returning NaNs.')
    stdResponse = NaN;
    steTrace = NaN;
end
