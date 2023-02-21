function [shore, bar, deep] = smoothInputs(shore,bar,deep, dy)
%   [shore, bar, deep] = smoothInputs(shore,bar,deep, dy);
%
% interpolate input data to a short spacing dy using a spline

shore.ys = [min(shore.y): dy: max(shore.y)];
shore.xs = interp1(shore.y,shore.x,shore.ys,'spline');

bar.ys = [min(bar.y): dy: max(bar.y)];
bar.xs = interp1(bar.y,bar.x,bar.ys,'spline');

deep.ys = [min(deep.y): dy: max(deep.y)];
deep.xs = interp1(deep.y,deep.x,deep.ys,'spline');
