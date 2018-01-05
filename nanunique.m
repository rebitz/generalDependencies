function [uniques] = nanunique(vect)

[uniques] = unique(vect);

uniques = uniques(~isnan(uniques));