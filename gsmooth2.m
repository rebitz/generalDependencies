function y=gsmooth2(x,sig,dimension)
% y=gsmooth(x,sig,dimension)
%   smooths row vector x by a gaussian filter with standard deviation of sig
%   
% if x is a matrix, smooth by rows, unless otherwise specified by dimension
%   dimension = 1: smooth by rows
%   dimension = 2: smooth by columns
%
% written long ago by JM Pearson
% adapted for use on matrixes by RB Ebitz 9/13

% try,
    
if nargin < 3
    dimension = 1; % default to smoothing across rows
end
%%     
if (size(x,1)==1) || (size(x,2)==1) % if we get a vector
    if(size(x,1)~=1) % but it's the wrong orientation
        x=x';
    end
    mat = 0;
    rows = 1;
else % if we get a matrix
    if dimension == 2 % but we want to smooth by colums
        x = x';
    end
    mat = 1;
    rows = size(x,1);
end
    
y = NaN(size(x,1),size(x,2));

for j = 1:rows
    
    xlocal = x(j,:);
    
    tot_length=size(xlocal,2); %total filter length

    filter=normpdf(-tot_length:1:tot_length,0,sig);

    ylocal=[];

    notnan=~isnan(xlocal);

    for i=1:tot_length
        index_range=-(i-1):1:(tot_length-i);
        index_range_shifted=index_range+tot_length+1;
        new_filter=filter(index_range_shifted); %translate gaussian and cut off at endpts
        new_filter=new_filter(notnan); % remove bins that will hit NaNs
        new_filter=new_filter/sum(new_filter); %normalize
        ylocal(:,i)=new_filter*xlocal(notnan)';
    end
    
    y(j,:) = ylocal;
end

if mat && dimension == 2;
    y = y';
end

% catch,
%     q=lasterror;
%     keyboard
% end