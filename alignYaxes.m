function alignYaxes(figHandle,ylimOpt)

if nargin < 1
    fH = gcf;
else
    fH = figHandle;
end

% subplotHandles
sH = get(fH,'Children');

if nargin < 2
    tmp = [inf -inf];
    for i = 1:length(sH)
        tmpLim = get(sH(i),'YLim');
        finalLim(1) = min([tmp(1), tmpLim(1)]);
        finalLim(2) = max([tmp(2), tmpLim(2)]);
        tmp = tmpLim;
    end
else
    finalLim = ylimOpt;
end

for i = 1:length(sH)
    set(sH(i),'YLim',[finalLim(1) finalLim(2)]);
end