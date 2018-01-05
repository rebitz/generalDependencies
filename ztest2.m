function [h,p,zStat] = ztest2(P1,P2,varargin)
%ZTEST2 two tailed, two proportion z-test, pooled variance estimate.
%   H = ZTEST2(P,N) performs a t-test of the hypothesis that two
%   independent samples of draws from binomial processes, in the data
%   vectors P1 and P2, come from distributions with equal probability of
%   sucess. H is the result of this statistical test, evaluated at p <
%   0.05.
%
%   P1 and P2 should be vectors of [proportion sucess, total observations]
%
%   The proportion can be expressed either as a count (i.e. 3 sucesses) or
%   a proportion (i.e. 0.3). The script treats input values > 1 as counts.
%
%   example useage: [h,p] = ztest2([0.1 100],[0.5 100])
%
%   Evaluates the hypothesis that a sucess rate of 0.1, observed on 100
%   independant draws, is significantly different from 0.5.
%
%   See also ZTEST, TTEST2.


if nargin < 2
    fprintf('\n');
    fprintf('Error: too few input arguments \n');
end

% Process remaining arguments
alpha = 0.05; % default to p < 0.05
tail = 0;    % code for two-sided;

% initialize
h = NaN; zStat = NaN; p = NaN;

if nargin>=3
    alpha = varargin{1};
end

if length(P1) == 2 && length(P2) == 2
    n = [P1(2),P2(2)];
    prop = [P1(1),P2(1)];
    
    % rescale to proportion if input seems to be a count
    if prop(1) > 1 || prop(2) > 1
        fprintf('\n');
        fprintf('Treating input proportions as counts. \n');
        prop  = prop ./ n;
    end
    
    if n(1) < 10 || n(2) < 10
        fprintf('\n');
        fprintf('Error: n observations must be greater than 10 \n');
        return
    end
    
    % 2 proportions z test:
    pooledP = (prop(1)*n(1)+prop(2)*n(2))./(n(1)+n(2));
    pooledSE = sqrt(pooledP*(1-pooledP)*((1./n(1))+(1./n(2))));
    zStat = (prop(1)-prop(2))/pooledSE;

else
        fprintf('\n');
        fprintf('Error: Wrong input format \n');
        fprintf('P1 and P2 must be 2x1 or 1x2 vectors, proportion and n \n');
        return
end

if zStat > 0
    p = (1-normcdf(zStat))*2;
else
    p = (normcdf(zStat))*2;
end

% high and low threshold values
thresh = [norminv(alpha/2,0,1), -norminv(alpha/2,0,1)];

if zStat < thresh(1) || zStat > thresh(2)
    h = 1;
else
    h = 0;
end