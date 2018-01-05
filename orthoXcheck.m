function [quality] = orthoXcheck(X)
%
% [quality] = orthoXcheck(predictorMatrix)
% 
% orthoXcheck implements Gram-Schmidt orthogonalization
%
% orthoXcheck takes a predictor matrix for a regression and checks to make
% sure there are no dependencies amongst the columns
%
% The method is taken from Applied Regression Analysis (3rd ed.), by
% Draper & Smith

%%

% % example from the text to make sure we're coding it correctly
% X = ones(5,1);
% X = [X [1:5]'];
% X = [X [1:5].^2'];
% X = [X [1:5].^3'];

% check for intercept in X mx, append if absent
if sum(X(:,1) == 1) ~= length(X(:,1))
    X = [ones(size(X,1),1),X];
end

% remove any rows w/ NaNs (unfortunately)
X = X(sum(isnan(X),2)~=0,:);

Z = X(:,1); % setup transformed mx

for i = 2:size(X,2) % for each column
    % b/c I fail at linear algebra:
    % M.^-1 = pinv(M)
    % M' = M'
    
    Z(:,i) = X(:,i) - Z*pinv(Z'*Z)*Z'*X(:,i);
end

% at this point Z is a matrix of orthogonally transformed predictors
% if any of these columns contains mostly 1s or v. small numbers, that
% indicates a dependency in the matrix where some of the predictors are
% linear combinations of each other

sum(Z==0)
