function [mu,idx] = bootstrpMean(nBoot,dat)

%%
idx = randi(length(dat),length(dat),nBoot);
mu = mean(dat(idx))';