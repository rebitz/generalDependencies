function [theta,Zlabels,lik] = exp3mixSpecial(times,rewards,theta)
% fit a mixture of 3 exponentials
%  theta is fitted parameters:
%   theta(1) hazard of fixed distribution
%   theta(2) hazard of rwd == 1
%   theta(3) hazard of rwd == 0
%   theta(4) mixing proportion

%%
epsilon = 10^-6;

% seed according to my hypothesis
t1 = nanmean(times(times >= (1/3)*max(times)+1));
t2 = nanmean(times(times <= (2/3)*max(times)+1));
t3 = nanmean(times(times <= (2/3)*max(times)+1));

% seed according to the opposite H1
% t1 = nanmean(times(times <= (2/3)*max(times)+1));
% t2 = nanmean(times(times >= (1/3)*max(times)+1));
% t3 = nanmean(times(times >= (1/3)*max(times)+1));

t4 = 0.5; % start by assuming an equal mixture
tHat = [t1, t2, t3, t4];

% setup our formula for the pdf
% D = @(x,rwd,tHat) (((tHat(1)/(1+tHat(1))).^x)*(tHat(4)/(1+tHat(1)))) + ...
%     (rwd .* (((tHat(2)/(1+tHat(2))).^x)*((1-tHat(4))/(1+tHat(2))))) + ...
%     ((1-rwd) .* (((tHat(3)/(1+tHat(3))).^x)*((1-tHat(4))/(1+tHat(3)))));
% 
% fZ1 = @(x,rwd,tHat) (((tHat(1)/(1+tHat(1))).^x)*(tHat(4)/(1+tHat(1)))) ./ D(x,rwd,tHat);
% fZ2 = @(x,rwd,tHat) (rwd .* (((tHat(2)/(1+tHat(2))).^x)*((1-tHat(4))/(1+tHat(2))))) ./ D(x,rwd,tHat);
% fZ3 = @(x,rwd,tHat) ((1-rwd) .* (((tHat(3)/(1+tHat(3))).^x)*((1-tHat(4))/(1+tHat(3))))) ./ D(x,rwd,tHat);

% setup our formula for the pdf
D = @(x,rwd,tHat) (((tHat(1)/(1+tHat(1))).^x)*(tHat(4)/(1+tHat(1)))) + ...
    (((tHat(2)/(1+tHat(2))).^x).*((rwd .* (1-tHat(4)))/(1+tHat(2)))) + ...
    (((tHat(3)/(1+tHat(3))).^x).*(((1-rwd) .* (1-tHat(4)))/(1+tHat(3))));

fZ1 = @(x,rwd,tHat) (((tHat(1)/(1+tHat(1))).^x)*(tHat(4)/(1+tHat(1)))) ./ D(x,rwd,tHat);
fZ2 = @(x,rwd,tHat) ((((tHat(2)/(1+tHat(2))).^x).*((rwd .* (1-tHat(4)))/(1+tHat(2))))) ./ D(x,rwd,tHat);
fZ3 = @(x,rwd,tHat) ((((tHat(3)/(1+tHat(3))).^x).*(((1-rwd) .* (1-tHat(4)))/(1+tHat(3))))) ./ D(x,rwd,tHat);


disp('starting EM');

i = 0; stop = false;
while ~stop

    % E-step
    Z1 = fZ1(times,rewards,tHat);
    Z2 = fZ2(times,rewards,tHat);
    Z3 = fZ3(times,rewards,tHat);

    % M-step
    t1 = sum(Z1.*times)/sum(Z1);
    t2 = sum(Z2.*times)/sum(Z2);
    t3 = sum(Z3.*times)/sum(Z3);
    t4 = nanmean(Z1);

    % update theta estimates for next round
    lastTHat = tHat;
    tHat = [t1 t2 t3 t4];
    
    % but stop if we've reached our stopping criteria
    stop = sum(and(1-epsilon <= tHat./lastTHat,tHat./lastTHat <= 1+epsilon)) == length(tHat);
    i = i+1;
end

fprintf('EM finished after %d iterations',i)
theta = tHat;
fprintf('\n mixture 1 hazard: %2.2f',theta(1))
fprintf('\n mixture 2 hazard: %2.2f',theta(2))
fprintf('\n mixture 3 hazard: %2.2f',theta(3))
fprintf('\n weights: %2.2f \n',theta(4))

% get the max lik label:
Z1 = fZ1(times,rewards,tHat);
Z2 = fZ2(times,rewards,tHat);
Z3 = fZ3(times,rewards,tHat);

[~,Zlabels] = max([Z1;Z2;Z3]);

% now the formula for the marginal
f3 = @(x,rwd) ((theta(4)).*(1./(1+theta(1))).*((theta(1)./(1+theta(1))).^x)) + ...
    (rwd .* (1-theta(4)).*(1./(1+theta(2))).*((theta(2)./(1+theta(2))).^x) + ...
    ((1-rwd) .* (1-theta(4)).*(1./(1+theta(3))).*((theta(3)./(1+theta(3))).^x))); % p 39 mT, first eq

lik = sum(log(f3(times,rewards)));

% % just skip this plot since it muddies the waters
% h = plot(x+1,f3(x),'-.k');
% set(h,'LineWidth',2)