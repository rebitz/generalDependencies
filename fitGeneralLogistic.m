function z = fitGeneralLogistic(params,t,y)
% params are starting parameters
% t = "time" points, or rather x data
% y = y data to be fit
% outputs z - which should be the cost of the function

r = params(1) % slope
t0 = params(2) % intercept
A = params(3) % minimum
K = params(4) % maximum

h = A + ((K-A)./(1+exp(-r*(t-t0))));

numerator = -dot(h,y)^2;
denominator = dot(h,h);
z = numerator/denominator

end