function [w,xpos] = wienerFilter(X,y,varargin)
%WIENERFILTER extract the best fit linear filter on X for
%   predicting y, controlling for autocorrelations in X
%
%   w = WIENERFILTER(X,y) returns the weights at each time lag
%
%   w = WIENERFILTER(X,y,nlags) allows you to specify the number of lags
%       default is 10 lags
%
%   [w,xpos] = WIENERFILTER(X,y,nlags) returns the x-position associated
%       with each of the weights
%
%   w = WIENERFILTER(X,y,[],'causal')
%       use the causal version of the filter? (one-sided, before the event)
%       else, the default version is the non-causal/symmetrical version
%
%   w = WIENERFILTER(X,y,[],'uncorrected')
%       don't correct for the autocorrelations in X


if nargin < 2
    error('not enough input arguments')
end

% set the defaults
nlags = 10;
applyCorrection = 1; % correct for autocorrelations in X?
causalFilter = 0; % only look at pre window,
% else, we'll work symmetrically

% change the defaults based on inputs
if length(varargin)>=1 % if the number of lags is present
    nlags = varargin{1};
end

if length(varargin)>1 % if we have more than 2 var args in    
    % find the pre-event, causal filter only
    causalFilter = sum(strcmp(varargin,'causal'))>0;
    
    % correct for autocorrelations, unless told otherwise
    applyCorrection = sum(strcmp(varargin,'uncorrected'))==0;
    
    % warnings for typos - if we got to here, but changed nothing
    if applyCorrection == 1 && causalFilter == 0
        disp('CAUTION: str inputs detected, but nothing was implemented');
        disp('  running default behavior: non-causal filter with a correlation correction');
        disp('  try again using the flags "causal" and/or "uncorrected"?');
    end
end

%% main body of the code

% transform input into col vectors if necessary
if length(X) > size(X,1)
    X = X';
end
if length(y) > size(y,1)
    y = y';
end

% make some temporary copies
Xtmp = X;
ytmp = y;

% center the vars
x = Xtmp - mean(Xtmp);
y = ytmp - mean(ytmp);

if causalFilter
    % autocorrelation
    [xc,lgs1] = xcorr(x,x,nlags,'biased');
    r = [xc(nlags+1:end)]; % only hang onto the pre-event
    rxx = toeplitz(r);

    % cross correlation
    [c,lgs2] = xcorr(x,y,nlags,'biased');
    p = c(1:nlags+1); % only hang onto the pre-event

    if applyCorrection
        w = pinv(rxx)*p;
    else
        w = p; % no controlling for autocorrelations
    end

    xpos = [-nlags:0]';
else
    % autocorrelation
    [xc,lgs1] = xcorr(x,x,nlags*2,'biased');
    r = [xc(nlags*2+1:end)];
    rxx = toeplitz(r); % make the toeplitz matrix

    % cross correlation
    [c,lgs2] = xcorr(x,y,nlags,'biased');

    % apply the correction
    if applyCorrection
        w = pinv(rxx)*c; % symmetrical
    else
        w = c;
    end

    xpos = lgs2';

end
