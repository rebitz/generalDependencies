function [h,resid] = olsfit(x,y)
%OLSLINE  Orthogonal least squares regression line.
%
%   OLSLINE(x,y) takes the x-plotted and y-plotted features of data and
%   fits a line in the total least squares sense, using unconstrained
%   minimization. You'd use this in place of lsline when comparing two
%   dependant variables, where there is error in measurement along both 
%   the x and y axes.
%
%   [h] = OLSLINE(x,y) gives the handle to the line.
%   [h,resid] = OLSLINE(x,y) gives sum of squared error of fit.

% sum of squared perpendicular distances
R = @(beta) sum((abs(beta(1)*x-y+beta(2))./sqrt(beta(1).^2+1)).^2);

% initial guesses - can be anything
beta0 = -.1; %offset
beta1 = 0; %slope
x0 = [beta1 beta0];

% do minimixation
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
[B,resid] = fminunc(R,x0,options);

% plot it in the current axis
hold on;
h = refline(B);
hold off;