function [pos] = findCentroid(positionVector)
%findClustroid finds the geometric centroid of a bunch of points
%   ROWS are observations
%   COLUMNS are measures w/in each dimension [x,y,z,etc]

goods = sum(isnan(positionVector),2)==0;
positionVector = positionVector(goods,:);

pos = mean(positionVector);

% ptDists = NaN(length(positionVector),1);
% for i = 1:length(positionVector)
%     pt = positionVector(i,:);
%     pt = repmat(pt,length(positionVector),1);
% 
%     % calc distance
%     dists = 0;
%     for j = 1:size(positionVector,2)
%         dists = dists + ((pt(:,j)-positionVector(:,j)).^2);
%     end
%     dists = sqrt(dists);
%     ptDists(i) = nansum(dists); % THIS IS WHERE YOU SPECIFY THE FUNCTION
%     % so currently finds the point that minimizes sum sq err
% end
% 
% [~,id] = min(ptDists);
% pos = positionVector(id,:);
            
end

