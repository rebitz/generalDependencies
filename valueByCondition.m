function fH = valueByCondition(rewards,choices,condition)
% takes cell arrays rewards and choices
%   plots p(choose|value)
% if given condition information, does this separtely for each condition
%   for 3 target case, uses the mean of the alternatives

%%
nBins = 10;

rwdDiff = [-1.1,1.1];
rwdBins = min(rwdDiff):range(rwdDiff)/nBins:max(rwdDiff);
% rwdBins(end) = rwdBins(end)+1;

prefCurves = NaN(2,nBins,length(rewards));

for k = 1:length(rewards);
    
    if ~isempty(condition{k})
        
    states = condition{k};
    rwds = rewards{k}./100;
    chs = choices{k};
    
    col = 1:length(chs);
    chosenValue = rwds(sub2ind(size(rwds),chs,col));
    
    % advance and modulo divide the index to get the other values
    otherValues = NaN(size(rwds)-[1 0]);
    for i = 1:size(rwds,1)-1
        plusOne = mod((chs+i), size(rwds,1));
        plusOne(plusOne==0) = size(rwds,1);
        otherValues(i,:) = rwds(sub2ind(size(rwds),plusOne,col));
    end
    
    % now, we want to re-express each value as - mean of the alternatives
    allValues = [chosenValue; otherValues];
    rwdDiff = NaN(size(allValues));
    
    for i = 1:size(rwds,1)
        tmp = 1:size(rwds,1);
        others = pop(tmp,i);
        rwdDiff(i,:) = allValues(i,:) - nanmean(allValues(others,:),1);
    end

    chs = [ones(size(chosenValue)); zeros(size(otherValues))];
    states = [repmat(states,size(rwds,1),1)];
    
    % flatten
    states = states(:); chs = chs(:); rwdDiff = rwdDiff(:);
    
    for i = 1:nBins
        binIdx = and(rwdDiff >= rwdBins(i), rwdDiff < rwdBins(i+1));
        
        % THIS IS FUCKED - GO BACK TO THE RECORDING STUFF!
        tmp = unique(states);
        for state = 1:length(unique(states));
            tmpIdx = and(binIdx,states==tmp(state));
            prefCurves(state,i,k) = nanmean(chs(tmpIdx));
        end
    end
    
    end
    
%     keyboard();
end

% plotz:

x = rwdBins(1:end-1)+mean(diff(rwdBins))/2;
% x = x/100;

m = nanmean(prefCurves,3)';
e = nanstd(prefCurves,[],3)'./sqrt(size(prefCurves,3)-1);

fH = figure(99); hold on;
set(gca,'FontSize',16);

h = plot(x,m,'.','MarkerSize',25,'LineStyle','none');
set(h(1),'Color','k'); set(h(2),'Color','b')
sigline;

legend(flipud(h),'explore','exploit','Location','NorthWest');

h = errorbar([x',x'],m,e);
for i = 1:length(h); removeEBarEnds(h(i)); end
set(h(1),'Color','k'); set(h(2),'Color','b');
set(h,'LineStyle','none');

ylabel('p(chose T1)');
xlabel('p(reward | T1) - p(reward | T2)');
