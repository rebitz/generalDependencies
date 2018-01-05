function out = ste(data)

out = std(data)./sqrt(sum(~isnan(data)));