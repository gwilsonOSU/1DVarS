function [xy, d] = findClosestPoint(a, xy0)
%   [xy, d] = findClosestPoint(abcLine, xy0)
%
%   find the point, xy, on a line that is closest to a query point, xy0,
%   and the orthogonal distance, d.  The line is defined as
%   ax+by+c=0.  abcLine is 1x3 and xy0 is 1x2.
%
%  equations are from wikipedia "distance from a point to a line"

ab2 = a(1)^2+a(2)^2;
d = abs(a*[xy0 1]')/sqrt(ab2);
foo = a(2)*xy0(1) - a(1)*xy0(2);
xy = [(a(2)*foo-a(1)*a(3)) (-a(1)*foo-a(2)*a(3))] / ab2;
