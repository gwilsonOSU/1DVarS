function h = makeBackgroundProfileComposite(x, xOff, hOff, betaShore, betaOff)
%
%   h = makeBackgroundProfileComposite(x, xOff, hOff, betaShore, betaOff)
%
%  Create a background beach profile of the form
%       h = a + beta*x + gamma*exp(-kx)
%  that merges between a planar offshore profile and a concave up nearshore
%  profile.  Parameters are constrained by an offshore depth and slope and
%  have a e-folding scale, k, that is related to betaShore, the foreshore
%  beach slope (Ruessink et al, 2003).  
%  This function is designed to provide a plausible underlying mean beach
%  profile to which sand bars can be added using the Ruessink algorithms
%  and the matlab file makeParametricBeach.
%  Inputs:
%   x           - regularly spaced cross-shore positions for estimation
%   xOff, zOff  - a match point offshore of the exponential
%   betaOff     - linear slope at the match point
%   betaShore   - slope of shoreshore

%  Holman, 09/11

% check if offshore slope intersects beach above z=0.  If so, settle with
% planar profile based on offshore knownn point.
if (betaOff*xOff > hOff)    % too steep offshore
    h = hOff/xOff * x;
    disp('WARNING - betaOff too steep.  Using plane beach approximation')
else
    [k, gamma] = findKGammaComposite(xOff, hOff, betaShore, betaOff);
    h = gamma * (exp(-k*x) - 1) + betaOff*x;
end

