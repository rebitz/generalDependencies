function y = generalLogistic(params,t)

r = params(1); % slope
t0 = params(2); % intercept
A = params(3); % minimum
K = params(4); % maximum

y = A + ((K-A)./(1+exp(-r*(t-t0))));

end