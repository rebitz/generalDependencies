function [p] = derrange(n)

p = 1:n;
for k = 1:n
   r = ceil(rand*n);
   while p(k) == r || p(r) == k % Reject non-derangement
      r = ceil(rand * n);
   end
   t = p(k);
   p(k) = p(r);
   p(r) = t;
end