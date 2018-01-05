function y = logistic(params,xdata)

y = 1./(1+exp(-params(1)*(xdata-params(2))));

end