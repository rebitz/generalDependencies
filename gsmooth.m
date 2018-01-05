function y=gsmooth(x,sig)
% y=gsmooth(x,sig)
%   smooths row vector x by a gaussian filter with standard deviation of sig
% 
%
% written long ago by JM Pearson

try,
    
    
if(size(x,1)~=1) % but it's the wrong orientation
    x=x';
end

tot_length=size(x,2); %total filter length

filter=normpdf(-tot_length:1:tot_length,0,sig);

y=[];

notnan=~isnan(x);

for i=1:tot_length
    index_range=-(i-1):1:(tot_length-i);
    index_range_shifted=index_range+tot_length+1;
    new_filter=filter(index_range_shifted); %translate gaussian and cut off at endpts
    new_filter=new_filter(notnan); % remove bins that will hit NaNs
    new_filter=new_filter/sum(new_filter); %normalize
    y(:,i)=new_filter*x(notnan)';
end

catch,
    q=lasterror;
    keyboard
end