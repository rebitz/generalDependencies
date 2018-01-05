function map = Heatmap(x, y, xlim, ylim, nbinx, nbiny)
% Heatmap
% map = Heatmap(x, y, xlim, ylim, nbinx, nbiny)
%     x and y are vectors of x,y coordinate pairs.
%     xlim and ylim define the boundaries of interest; values outside these
%         ranges will be discarded. Should have the form [xmin xmax].
%     nbinx and nbiny are the number of bins.
%     map is a matrix of size [nbiny nbinx] which counts the number of
%         occurrences in each bin (ie, a 2-d histogram)
% 
% Heatmap(data(i).Reyepts(1),data(i).Reyepts(2),xedges,yedges,100,100);
% 

global data i

xmin = xlim(1);
xmax = xlim(2);
ymin = ylim(1);
ymax = ylim(2);

% Identify which points have x OR y values that violate either boundary:
outOfBounds = x<=xmin | x>xmax | y<=ymin | y>ymax;
% Throw away the ones that do:
x = x(~outOfBounds);
y = y(~outOfBounds);

x = x - xmin;
y = y - ymin;
% Now x ranges between 0 and xmax-xmin, likewise for y

% Remap x and y onto [0 nbinx] and [0 nbiny]:
x = x * nbinx / (xmax - xmin);
y = y * nbiny / (ymax - ymin);

% x and y are continuous values. We need to convert them into grid indices.
x = ceil(x);
y = ceil(y);

% Now, we have row, column indices. What is the easiest way to count how
% many repetitions we have? We can convert our row, column indices into
% vector indices by using sub2ind.
ind = sub2ind([nbiny nbinx], y, x);

% Count the number of times each ind occurs:
n = histc(ind,1:nbiny*nbinx);

% Produce the heat map:
map = reshape(n, nbiny, nbinx);