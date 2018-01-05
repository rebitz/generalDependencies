function [M,E,xpos,count] = extractErrorbar(xdata,ydata,nbins,dist)
%
% [M,E,xpos] = extractErrorbar(xdata,ydata,nbins,distribution)
%
% breaks up ydata by dividing associated xdata into bins
% 
% returns yhat for that bin as M (mean), E (ste), and xpos (center of bins)
%   could also return raw count [M,E,xpos,count]
%
% nbins = number of bins to divide into
%   nbins can be set high to get only unique values of x
%
% distribution = controls type of error bar calculation
%   if 'normal' does standard normal distribution errors (std/sqrt(sum-1))
%       (default is 'normal')
%   if 'binomial' does binomial standard errors (sqrt(p*(1-p))/sum)
%

if nargin < 3 || isempty(nbins)
    nbins = 10;
end

if nargin < 4
    dist = 'normal';
end

%%

if length(unique(xdata)) > nbins
    bat = quantile(xdata,[0:1/nbins:1]);
    bat(1) = min(xdata)-1;

    [M,E,count] = deal(NaN(length(bat)-1,1));

    for i = 1:length(bat)-1

        idx = and(xdata > bat(i), xdata <= bat(i+1));

        M(i) = nanmean(ydata(idx));
        count(i) = sum(idx);
        
        if strcmp(dist,'normal')  
            E(i) = nanstd(ydata(idx)) ./ sqrt(sum(idx)-1);
        elseif strcmp(dist,'binomial');
    %         [~,V] = binostat(sum(idx),M(i))
            E(i) = sqrt(M(i)*(1-M(i))/sum(idx));
        end
    end

    xpos = bat(1:end-1)+diff(bat)/2;
else
    disp('too few data points for quantiles!')
    bat = unique(xdata);

    [M,E,count] = deal(NaN(length(bat),1));

    for i = 1:length(bat)

        idx = xdata == bat(i);

        M(i) = nanmean(ydata(idx));
        count(i) = sum(idx);

        if strcmp(dist,'normal')  
            E(i) = nanstd(ydata(idx)) ./ sqrt(sum(idx)-1);
        elseif strcmp(dist,'binomial');
    %         [~,V] = binostat(sum(idx),M(i))
            E(i) = sqrt(M(i)*(1-M(i))/sum(idx));
        end
    end

    xpos = bat;
end