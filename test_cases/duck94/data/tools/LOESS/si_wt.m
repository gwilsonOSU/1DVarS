function w = si_wt(x, w)
%
% [a] = si_wt(x, w)
%
% Input
%   x are nxm inputs, weights will be centered on x=0
%   w are nx1 values of Root normalized signal/total variance of data
%
% Output
%   a are weights 0<=w<=1

[n,m] = size(x);
LN = 2; % fix up so scaling sends weights to zero at x=1
x = x*LN;
if(~exist('w'))
    w = ones(n,1);
else
    w=(w+1e-3)/(1+1e-3);
end
% compute data-data covariance matrix (true variance at data points + noise along% diagonal)
% assume true variance is constant at all locations (var=1)
r = 0;
for i=1:m
    r = r + (repmat(x(:,i)', n, 1) - repmat(x(:,i), 1, n)).^2;
end
r = exp(-r).*(w*w'); % has covariance of w.^2 along diagonal
% diagonal elements have observation error added
%r = r + diag(1-w); % not quite right, if w=0, get cov of 2 on diagnal. should be one, and zero everywhere else
r = r + diag(1-w.^2); % now, if r is scaled correctly, this oughta work

% invert it
P = ceil(95 + 5*(100/(100+n)));
%P=100;% only reduce P if many data

r = svd_invP(r, P); % approach 95% for large N
%r = inv(r);

% now, data-true covariance vector
ri = exp(-(x.^2)*ones(m,1)).*(w.^2);

% calculate the weights
w = r*ri;