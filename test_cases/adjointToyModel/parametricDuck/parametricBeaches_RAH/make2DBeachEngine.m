function [x,y,h] = make2DBeachEngine(shore,bar,deep, hSea, x, y)
%   [x,y,h] = make2DBeachEngine(shore,bar,deep,hSea,{xGrid yGrid})
%
%   Create a synthetic 2DH bathymetry using a 2D implementation of the 1D
%   parametric barred beach profile algorithm detailed in Holman et al,
%   2014, Coastal Engineering.
%   At present, x is expected to be roughly offshore and y alongshore.  
%   The user supplies approximate data structures to describe the shoreline,
%   a sand bar and a deep water depth contour as follows:
%  shore
%   - x,y, a set of x-y points representing the approx shoreline (defined by
%       a user-selected datum.  This can be coarse since it will be
%       smoothed and interpolated in the subsequent processing
%   - yb, beta, the climatological beach slope over the same shoreline extent,
%       specified points yb.  These slopes should be a guess at the long
%       term average slope and will often be a constant value, but
%       represented at a set of yb points that might be the same as the y
%       points above (but must be specified).
%
%  bar
%   - is a structure with only a reasonable set of x,y locations that
%       correspond to the crest of a sand bar.  This need not be the inner
%       (first) bar, just one that is convenient (perhaps shows up easily
%       in remote sensing)
%
%  deep
%   - a structure that specifies any deep water line of bathymetry that is
%       outside (deeper than) the zone of wave-generated sand bars.  This
%       might often be a 15-20 m depth contour.  
%       The structure has fields x,y that can be specified coarsely (and
%       will be subsequently interpolated and smoothed.  At each of these
%       points the user also must specify a depth, h, and approximate
%       slope, beta.
%   hSea, max depth of bar activity, empirical, from Ruessink paper.
%
%   The user may optionally specify the x and y values of a regular grid at
%   which bathymetry will be estimated.  If not specified, the algorithm
%   will return a sensible set of axes.
%
%   The routine returns the x and y grid axes as well as bathymetry, h, as
%   a 2D matrix of size length(y) by length(x).

% start by improving the inputs by shorting them into an increasing order,
% then smoothing and interpolating.  Results are stored in the same
% structures but with fields xs,ys for smoothed values

[shore, bar, deep] = sortInputs(shore,bar,deep); % ensure increasing order
% smooth the shore, bar and deep to short dy spacing
dy = 1;             % 1 m spacing is default but this may be over fine
[shore, bar, deep] = smoothInputs(shore,bar,deep,dy);
    
% define a map grid.  If the user supplies this, fine, otherwise supply a
% reasonable alternative.  Note that the user grid can include dry beach
% regions
if (nargin == 3)    % no grid supplied
    [x,y] = designGrid(shore,deep);
end
[X,Y] = meshgrid(x,y);
X = X(:); Y = Y(:);
xMin = 4;           % only analyze if at least this far offshore.
h = nan(size(X));

% now consider each point in turn.
for i = 1: length(X)        % loop through all points
    handle = waitbar(i/length(X));
    % for each point, find the mean shoreline over a span based on the
    % x-distance along with cross-shore transect and support result
    if (Y(i) < shore.ys(1)) || (Y(i) > shore.ys(end)) % out of y-range
        h(i) = nan;
    elseif ((X(i)-interp1([shore.ys], [shore.xs], Y(i))) < xMin)  % if on land, bail out
        h(i) = nan;
    else
        trans = findShoreNormalSmooth(X(i),Y(i),shore);
        trans = findBarAndDeep(trans, bar, deep);
        h(i) = findDepthOnTransect(trans, hSea);
    end
end

% finally reshape to a grid and delete the waitbar.
h = reshape(h,length(y), length(x));
delete(handle)
