function [sparsity] = giniSparsity(c)
%giniSparsity(c) calculate the Gini sparsity of vector c

N = length(c);

[~,sortIndx] = sort(c);

% calculate percent of total contributed by each:
pcnt = (c(sortIndx) ./ norm(c,1));

normX = ((N - [1:N] + (1/2))) ./ N;

sparsity = 1-2*sum(pcnt.*normX);

end

