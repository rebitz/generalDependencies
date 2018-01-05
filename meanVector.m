function [m,e,xpos,p,variance] = meanVector(y,x)
% finds the average and error of y, given some labelling scheme x

xpos = nanunique(x);

[m,e] = deal(NaN(1,length(xpos)));
for i = 1:length(xpos);
    m(i) = nanmean(y(x == xpos(i)));
    e(i) = nanstd(y(x == xpos(i)))./sqrt(sum(x == xpos(i))-1);
    variance(i) = nanvar(y(x == xpos(i)));
end

if length(xpos)==2
    [h,p] = ttest2(y(x == xpos(1)),y(x == xpos(2)))
else
    p = NaN;
end

end

