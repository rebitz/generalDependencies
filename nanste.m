function y = nanste(in)
%NANSTE Standard Error, ignoring NaNs

y = nanstd(in)./sqrt(sum(~isnan(in))-1);
