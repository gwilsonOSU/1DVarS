function r = hanning_wt(r)% % function r = hanning_wt(r)% % w = hanning_wt(x)% Input%   x are nxm inputs, weights will be centered on x=0%% Output%   w are weights 0<=w<=1% [n,m] = size(r);% convert to radial distancer = sqrt((r.^2)*ones(m,1));% hanning windowr = (1 - cos( pi*(0.5 + 0.5*r) ).^2).*(r<=1);