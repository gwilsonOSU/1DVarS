function t = findShoreNormalSmooth(x,y,s)
% trans = findShoreNormalSmooth(x,y,shore)
%
% given a structure, shore, with fields x,y,beta, yb, and an offshore point
% x,y, find the equation of a line y = mx + b that is orthogonal to the
% least squares shoreline fit over y-range that is +/-K times the mean
% offshore distance (done iteratively)

% returns transect structure, trans, with fields
%   trans.m     - slope of x-shore transect, i.e. y = mx + b
%   trans.b     - intercept of same
%   trans.xy    - closest point of mean shoreline to offshore query point
%   trans.d     - distance from query point to closest shoreline point
%   trans.b0    - interpolated foreshore slope at shore point xy.

K = 1.0;        % reasonable default value

nMin = 2;   % converged if this close to same number of points
dx = x - interp1([s.ys],[s.xs],y); 
oldInd = 1:length(s.ys);
ind = find((s.ys>y-K*dx) & (s.ys<y+K*dx));
while (abs(length(ind)-length(oldInd))>nMin)
    oldInd = ind;
    dx = x - mean(s.xs(ind));
    ind = find((s.ys>y-K*dx) & (s.ys<y+K*dx));
end

% fit shoreline using x locations as a function of y, then transform to y =
% mx+b format.  Note that this is more stable that direct solution of
% latter.
p = polyfit(s.ys(ind),s.xs(ind),1);
m = 1/p(1);
if isinf(m)
   m = 10^3;         % 89.94 degrees, close enough.
end
t.m = -1/m;                    % slope, int for xshore transects
t.b = y - t.m.*x;
[t.xy0,t.d] = findClosestPoint([1/p(1) -1 -p(2)/p(1)], [x y]);

% interp slopes, beta
if (t.xy0(2) < s.yb(1))
    t.beta0 = s.beta(1);
elseif (t.xy0(2) > s.yb(end))
    t.beta0 = s.beta(end);
else
    t.beta0 = interp1([s.yb],[s.beta],t.xy0(2)); % shoreline slopes
end
