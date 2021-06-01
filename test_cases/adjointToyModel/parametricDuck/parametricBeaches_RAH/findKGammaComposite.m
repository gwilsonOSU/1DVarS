function [k, gamma] = findKGammaComposite(xOff, hOff, betaShore, betaOff)
%
%   [k, gamma] = findKGammaComposite(xOff, hOff, betaShore, betaOff)
%
%  Routine to solve for the two unknowns of the exponential profile fit
%  based on data constraints.  Results and inputs are tied to the
%  parametric beach effort with sand bar profiles suggested by Ruessink et
%  al, 2003 and the results are usually input to the function
%  makeBackgroundProfileExp.  
%  Inputs:
%   xOff, hOff  - x,h values for some known point on the mean profile
%                 (usually beyond the active bar depth)
%   betaShore   - shoreline beach slope
%   beta        - deep water (outside bar zone) beach slope
%  Outputs:
%   k           - e-folding scale of nearshore exponential part of profile
%   gamma       - coefficient in exponential model.

% Holman, 10/11

% We require a numerical solution, done in the while loop below.

k = 1/xOff;         % first guess at k
fact = 0.1;
tol = 0.00001;                   % tolerance on search (in meters depth)
dh2 = (betaOff-betaShore)/k * (exp(-k*xOff)-1) + betaOff*xOff - hOff;
dh = (betaOff-betaShore)/(k*(1+fact)) * (exp(-(k*(1+fact))*xOff)-1) + betaOff*xOff-hOff;
if abs(dh)>abs(dh2)   % oops, wrong sign of search
    fact = -fact;
end

while (abs(dh)>tol)
    k = k*(1+fact);
    dh = (betaOff-betaShore)/k * (exp(-k*xOff)-1) + betaOff*xOff-hOff;
    if abs(dh)>abs(dh2)
        fact = -fact*0.1;
    end
    dh2 = dh;
end
gamma = (betaOff-betaShore)/k;
