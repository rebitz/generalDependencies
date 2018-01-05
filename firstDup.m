function vdup = firstDup(M)
    [SM Ord] = sort(M, 2);    % sort by row
    [rows cols] = find(~diff(SM, 1, 2));   % diff each row, and find indices of the repeated elements in sorted rows
    Mask = (size(M,2) + 1) * ones(size(M)); % create a Mask matrix with all size(M,2)+1
    ind = sub2ind(size(Ord), rows, cols+1); % add 1 to the column indices
    Mask(ind) = Ord(ind);   % get the original indices of each repeated elements in each row
    vdup = min(Mask, [], 2); % get the minimum indices of each row, which is the indices of first repeated element