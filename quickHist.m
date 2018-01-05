function h = quickHist(inputVector,currentAx,nbins)
%quickHist generate a quick histogram plot of a vector
%   Makes a plot 

if nargin < 2
    currentAx = 0; % default to new figure win
end

if nargin < 3
    nbins = 20;
end


tmp = inputVector;
xpos = [min(tmp):range(tmp)/nbins:max(tmp)];
if ~currentAx; figure(); end
xpos = xpos(1:end-1)+(diff(xpos)./2); % center the bins
h = bar(xpos,histc(tmp,xpos),'FaceColor','k');

end

