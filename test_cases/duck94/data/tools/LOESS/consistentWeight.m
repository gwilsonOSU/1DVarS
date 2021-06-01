function [w, s2] = consistentWeight(z, e2, wtol)
% compute weights consistent with observations and expected error variance for each observation
% 
% [w, s2] = consistentWeight(z, e2, wtol)
%
% Input
%   z, the observations
%   e2, the error VARIANCE (SQUARED, DAMMIT!) at each observation
%   wtol, (optional) convergence criteria (suggest 0.1 to 0.01)
%
% Output
%  w, the consistent weights
%  s2, the estimated true variance
%
% assume (wj^2) = s2/(s2 + e2j), where s2 is true variance and e2j is error variance at observation j
% solve problem by guessing wj, then compute s2, then update wj
% assumes that true variance, s2, is constant over the data

% if no wtol, specify
if (nargin < 3) | (nargin == 3 & isempty(wtol)) %(~exist('wtol')==1)
    tol = 0.01;
end

% initial guess, weights are uniform
N = length(e2);
winit = ones(N,1);
cnt = 0;
w = 2; % start to enter while loop
while (cnt <10 & max(abs(w-winit))>wtol)
    cnt = cnt + 1;
    % don't let large weights dominate too soon
    if(cnt>1)
        winit = (w+winit)/2;
    end
    % compute terms
    % weighted mean
    mu = sum(winit.*z)/sum(winit);
    nu = length(z); 
    % weighted mean squared residual-> true variance if weigths are accurate
    s2 =sum( (winit.*(z-mu)).^2 )/nu;
    % when nu is small, this estimate ought to be replaced by the a priori error (at least)
    s2 = mean(((nu-1)*s2 + e2)/nu);
    w = sqrt(s2./(s2 + e2));
    % no further convergence if nu=1
    if(nu==1)
        return
    end
end
