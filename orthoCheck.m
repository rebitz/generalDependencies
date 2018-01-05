function [orthoX,keptIdx] = orthoCheck(X,silent)
%
% orthoCheck takes a predictor matrix for a regression and checks the
% linear dependencies amongst the columns
%
% [orthoX] = orthoCheck(predictorMatrix)
%
% by default, this spits lots of information about how orthogonal the
% predictors are, but when optional argument "silent" == 1, it does not, eg
%   [orthoX] = orthoCheck(predictorMatrix,1) % runs silently
% 
% orthoCheck(X) implements Gram-Schmidt orthogonalization
%  returns orthoX, which is the Gram-Schmidt orthogonalization of the
%  predictor matrix X
% CAUTION: this method is highly sensitive to the order of regressors
%
% The method is taken from Applied Regression Analysis (3rd ed.), by
% Draper & Smith
%
% rbe, 7/28/14

%%
if nargin < 2
    % defaults to lots of feedback
    silent = 0;
end

% % example from the text to make sure we're coding it correctly
% X = ones(5,1);
% X = [X [1:5]'];
% X = [X [1:5].^2'];
% X = [X [1:5].^3'];

% check for intercept in X mx, append if absent
if sum(X(:,1) == 1) ~= length(X(:,1))
    intercepted = 1;
    X = [ones(size(X,1),1),X];
else
    intercepted = 0;
end

% remove any rows w/ NaNs (unfortunately)
keptIdx = (sum(isnan(X),2)==0);
X = X(sum(isnan(X),2)==0,:);

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

if ~silent
    
    fprintf('percent zero elements in each column')
    sum(Z==0)./size(Z,1)

    fprintf('percent v. v. small elements (0 +/- 0.001) in each column')
    sum(abs(0-Z)<0.001)./size(Z,1)

    fprintf('percent v. small elements (0 +/- 0.01) in each column')
    sum(abs(0-Z)<0.01)./size(Z,1)

    fprintf('percent small elements (0 +/- 0.1) in each column')
    sum(abs(0-Z)<0.1)./size(Z,1)

    % small is qualitative, relative to (0,1)
end


% % trying for principle component decomposition:
% Z = X;
% % % remove any nan's
% % Z = Z(sum(isnan(Z),2)==0,:);
% % 
% % % skip the intercept, but
% % for i = 2:size(X,2) % column equilibrate X
% %     Z(:,i) = Z(:,i)./sqrt(nansum((Z(:,i)-nanmean(Z(:,i))).^2));
% %     % each column should have sumsqr unity
% % end
% 
% % singular value decomposition:
% % [U,D,V] = svd(Z);
% % Z = U*D;
% 
% % all of the above should be equivalent to:
% [coeff,Z] = pca(Z,'Algorithm','svd')
% % but it is definitely not equivalent
% 
% % keyboard();

% remove the intercept if we added one
if nargout > 0
    if ~intercepted
        orthoX = Z;
    else
        orthoX = Z(:,2:end);
    end
end


