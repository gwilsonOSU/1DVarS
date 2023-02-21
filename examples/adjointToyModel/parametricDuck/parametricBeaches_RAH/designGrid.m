function [x,y] = designGrid(shore,deep)
%   [x,y] = designGrid(shore,deep)
%
% select x and y vectors of locations for bathy estimation based on the
% span of data supplied in the shoreline and deep vectors of input data.

% design specs
dx = 10;         % grid spacing in x and y
dy = 20;
incx = 100;      % this gives a rounding buffer around the data
incy = 100;      % in x and y

% find limits and make vectors
x0 = floor(min(shore.x)/incx) * incx;   % go from min shore to max deep
x1 = ceil(max(deep.x)/incx) * incx;
x = [x0:dx:x1];

y0 = floor(min(shore.y)/incy) * incy;   % from min to max y
y1 = ceil(max(shore.y)/incy) * incy;
y = [y0:dy:y1];
