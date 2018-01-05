function [seq] = fibIT(digits)
%fibIT generates a Fibonacci sequence, up to n digits

if nargin == 0
    digits = 10;
end

seed = [0 1];

seq = [seed, NaN(1,digits-length(seed))];

for i = length(seed):digits
    seq(i+1) = seq(i)+seq(i-1);
end

end

