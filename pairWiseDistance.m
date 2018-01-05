function [distance] = pairWiseDistance(positionVector)
%pairWiseDistance finds the total distance between points
%   pairWiseDistance(positionVector)
%       positionVector is rows = points, columns = measurements

% select rows without missing data
goods = sum(isnan(positionVector),2)==0;
positionVector = positionVector(goods,:);

ptDists = NaN(length(positionVector),1);
for i = 1:length(positionVector)
    pt = positionVector(i,:);
    pt = repmat(pt,length(positionVector),1);

    % calc distance
    dists = 0;
    for j = 1:size(positionVector,2)
        dists = dists + ((pt(:,j)-positionVector(:,j)).^2);
    end
    dists = sqrt(dists);
    
    ptDists(i) = nansum(dists); % THIS IS WHERE YOU SPECIFY THE FUNCTION
end

% sum of squared distances
distance = sum(ptDists.^2);

end
