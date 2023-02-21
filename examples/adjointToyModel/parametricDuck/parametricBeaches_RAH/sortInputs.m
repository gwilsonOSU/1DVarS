function [s,b,d] = sortInputs(s,b,d)
%   [shore,bar,deep] = sortInputs(shore,bar,deep)
%
% take input structures for shore, bar and deep and ensure they are in
% order of increasing y.  This just simplifies later algorithms.

[s.y,ind] = sort(s.y);
s.x = s.x(ind);
[s.yb,ind] = sort(s.yb);
s.beta = s.beta(ind);

[b.y,ind] = sort(b.y);
b.x = b.x(ind);

[d.y,ind] = sort(d.y);
d.x = d.x(ind);
d.h = d.h(ind);
d.beta = d.beta(ind);