function [theta,Zlabels,lik] = exp4mixSpecial(times,rewards,theta)
% fit a mixture of 3 exponentials
%  theta is fitted parameters:
%   theta(1) hazard of long,rwd == 1
%   theta(2) hazard of long,rwd == 0
%   theta(3) hazard of short,rwd == 1
%   theta(4) hazard of short,rwd == 0
%   theta(5) mixing proportion

%%
epsilon = 10^-6;

t1 = nanmean(times(times >= (1/3)*max(times)+1));
t2 = nanmean(times(times >= (1/3)*max(times)+1));
t3 = nanmean(times(times <= (2/3)*max(times)+1));
t4 = nanmean(times(times <= (2/3)*max(times)+1));
t5 = 0.5; % start by assuming an equal mixture
tHat = [t1, t2, t3, t4, t5]

% setup our formula for the pdf
D = @(x,rwd,tHat) (rwd.* ((tHat(1)/(1+tHat(1))).^x)*(tHat(5)/(1+tHat(1)))) + ...
    ((1-rwd).* ((tHat(2)/(1+tHat(2))).^x)*(tHat(5)/(1+tHat(2)))) + ...
    (rwd .* (((tHat(3)/(1+tHat(3))).^x)*((1-tHat(5))/(1+tHat(3))))) + ...
    ((1-rwd) .* (((tHat(4)/(1+tHat(4))).^x)*((1-tHat(5))/(1+tHat(4)))));

fZ1 = @(x,rwd,tHat) (rwd .* (((tHat(1)/(1+tHat(1))).^x)*(tHat(5)/(1+tHat(1))))) ./ D(x,rwd,tHat);
fZ2 = @(x,rwd,tHat) ((1-rwd) .* (((tHat(2)/(1+tHat(2))).^x)*((1-tHat(5))/(1+tHat(2))))) ./ D(x,rwd,tHat);
fZ3 = @(x,rwd,tHat) (rwd .* (((tHat(3)/(1+tHat(3))).^x)*((1-tHat(5))/(1+tHat(3))))) ./ D(x,rwd,tHat);
fZ4 = @(x,rwd,tHat) ((1-rwd) .* (((tHat(4)/(1+tHat(4))).^x)*((1-tHat(5))/(1+tHat(4))))) ./ D(x,rwd,tHat);

disp('starting EM');

i = 0; stop = false; maxIter = 5000;
while ~stop && i < maxIter

    % E-step
    Z1 = fZ1(times,rewards,tHat);
    Z2 = fZ2(times,rewards,tHat);
    Z3 = fZ3(times,rewards,tHat);
    Z4 = fZ4(times,rewards,tHat);

    % M-step
    t1 = sum(Z1.*times)/sum(Z1);
    t2 = sum(Z2.*times)/sum(Z2);
    t3 = sum(Z3.*times)/sum(Z3);
    t4 = sum(Z4.*times)/sum(Z4);
    t5 = nanmean(sum([Z1;Z2]));

    % update theta estimates for next round
    lastTHat = tHat;
    tHat = [t1 t2 t3 t4 t5];
    
    % but stop if we've reached our stopping criteria
    stop = sum(and(1-epsilon <= tHat./lastTHat,tHat./lastTHat <= 1+epsilon)) == length(tHat);
    i = i+1;
end

fprintf('EM finished after %d iterations',i)
theta = tHat;
fprintf('\n mixture 1 hazard: %2.2f',theta(1))
fprintf('\n mixture 2 hazard: %2.2f',theta(2))
fprintf('\n mixture 3 hazard: %2.2f',theta(3))
fprintf('\n mixture 4 hazard: %2.2f',theta(4))
fprintf('\n weights: %2.2f \n',theta(5))

% get the max lik label:
Z1 = fZ1(times,rewards,tHat);
Z2 = fZ2(times,rewards,tHat);
Z3 = fZ3(times,rewards,tHat);
Z4 = fZ4(times,rewards,tHat);

[~,Zlabels] = max([Z1;Z2;Z3;Z4]);

% now the formula for the marginal
f4 = @(x,rwd) (rwd .* (theta(5)).*(1./(1+theta(1))).*((theta(1)./(1+theta(1))).^x)) + ...
    ((1-rwd) .* (theta(5)).*(1./(1+theta(2))).*((theta(2)./(1+theta(2))).^x)) + ...
    (rwd .* (1-theta(5)).*(1./(1+theta(3))).*((theta(3)./(1+theta(3))).^x) + ...
    ((1-rwd) .* (1-theta(5)).*(1./(1+theta(4))).*((theta(4)./(1+theta(4))).^x))); % p 39 mT, first eq

lik = sum(log(f4(times,rewards)));

% % just skip this plot since it muddies the waters
% h = plot(x+1,f3(x),'-.k');
% set(h,'LineWidth',2)