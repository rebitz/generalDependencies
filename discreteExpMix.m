function [theta,Zlabels,lik] = discreteExpMix(x,K)
% DISCRETEEXPMIX fit a mixture of discrete exponential distributions
%   THETA = discreteExpMix(X,K) fits a mixture of K discrete exponenetial
%   distributions to the inter-event intervals in X (default, K=2). Returns
%   THETA, a vector of fitted model parameters. Entries of THETA are:
%       theta(1) survival parameter of first mixture (beta_1)
%       theta(2) weight of first mixture (pi_1)
%       theta(3:4) same for second mixture (beta_2, pi_2)
%       theta(5:6) same for third mixture (beta_3, pi_3), etc.
%   
%   [THETA,Zlabels] = discreteExpMix(X,K) also returns the most probable
%   mixing component for each observation in X (Z labels for mixtures 1:K)
%
%   [THETA,Zlabels,loglikelihood] = discreteExpMix(X,K) returns the
%   marginal logliklihood of the observations, p(X|THETA)

% this method comes out of this masters thesis:
%   https://ecommons.cornell.edu/handle/1813/2953
% "Mixtures of Exponential Distributions to Describe the Distribution of
% Poisson Means in Estimating the Number of Unobserved Classes", Barger

if min(x)~=0
    x = x - 1;
    subtracted = true;
else
    subtracted = false;
end

plotIt = true;
verbose = true;

if nargin < 2
    K = 2; % number of mixing components
end

%%
epsilon = 10^-6;
maxiter = 10^4;

% preallocate for parameters
[betas,weights] = deal(NaN(1,K));
[tHat] = deal(NaN(1,K*2));

% first, find some starting values by splitting up the data
bins = min(x):(1/(K+1))*max(x):max(x);
for k = 1:K
    betas(k) = nanmean(x(and(x >= bins(k), x <= bins(k+2))));
    weights(k) = (1/K); % start by assuming an equal mixture
end

% fill into our parameter vector
tHat(1:2:end) = betas;
tHat(2:2:end) = weights;

% now the discrete exponential (geometric) pdf:
%   geometric pdf: p * (1-p)^x
%     where p = 1./(1+beta)
geoF = @(x,beta) (1./(1+beta)).*(beta./(1+beta)).^x;

% now a formula for our weighted distributions:
wF = @(x,theta) theta(2)*geoF(x,theta(1));

%%

if K ~= 1
    disp('starting EM');
    exitFlag = true; % switches to false if we exit for a bad reason

    i = 0; stop = false;
    while ~stop

        % E-step: calculate the expectation of our observations | params
        [Z,pdf] = zLabels(x,tHat);

        % M-step:
        betas = tHat(1:2:end);
        weights = tHat(1:2:end);
        [betasN,weightsN] = deal(NaN(size(betas)));

        % update the betas:
        for k = 1:K
            betasN(k) = sum(Z(k,:).*x)/sum(Z(k,:));
        end

        % finally the weights, these sum to 1, so df for weights = K-1
        weightsN = mean(Z,2)';

        % update theta estimates for next round
        lastTHat = tHat;
        tHat(1:2:end) = betasN;
        tHat(2:2:end) = weightsN;    

        % but stop if we've reached our stopping criteria
        stop = sum(and(1-epsilon <= tHat./lastTHat,tHat./lastTHat <= 1+epsilon)) == length(tHat);
        stop = or(stop,i > maxiter);

        if sum(isnan(tHat)) > 0; stop = 1; exitFlag = false; end
        if verbose && round(i/100)*100 == i; lik = sum(log(sum(pdf))); fprintf('\n iter: %d, lik: %2.4f',i,lik); end

        i = i+1;
    end
else % if K does == 1
    % then our ML estimate is just mean of observations:
    tHat(1) = mean(x);
    tHat(2) = 1;
    exitFlag = true; i = 0;
end

if exitFlag
    fprintf('\n EM finished after %d iterations',i);
    
    theta = tHat;

    % get max lik mixture
    [Z,pdf] = zLabels(x,theta);
    [~,Zlabels] = max(Z);

    % likelihood:
    lik = sum(log(sum(pdf,1)));

    if plotIt
        figure();
        counts = histc(x,[0:max(x)]);
        semilogy([0:max(x)],counts/sum(counts),'.k','MarkerSize',20); hold on;
        [~,pdf] = zLabels([0:max(x)],theta);
        plot([0:max(x)]+0.5,sum(pdf,1),'LineWidth',2)
    end

    % if we subtracted 1 to make the data start at 0 to fit the model
    % we need to correct the survival parameter for interpretability
    if subtracted
        theta(1:2:end) = theta(1:2:end) + 1;
    end

else
    fprintf('\n EM failed to converge! \n');
end


function [Z,pdf] = zLabels(x,theta)
    % generate the zlabels for this iteration
    Z = NaN(K,length(x)); idx = [1,2];
    for k = 1:K
        Z(k,:) = wF(x,theta(idx));
        idx = idx+2;
    end
    pdf = Z;
    G = sum(Z,1); % sum over rows of Z
%     G = G + 10.^-64; % avoid divide by 0's
    Z = Z./repmat(G,K,1);
    
    %  zlabels are of the form:
    %     z1 = (pi_1*f_1(x,theta_1))./sum_n((pi_n*f_n(x,theta_n)))
end

end

