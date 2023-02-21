function h = makeBarredBeach(hShore, hSea, x, h0, xb)
%   hBarred = makeBarredBeach(hShore, hSea, x, hBackground, xb)
%
%  Create a parametric beach profile [x,hBarred] given a mean profile 
%  [x, hBackground] and
%  depths of the shoreward and seaward ends of the barred profile envelope,
%  hShore and hSea and the cross-shore phase of the bar.  xb is a bar crest
%  location which sets the phase, psi.  
%  This routine estimates a reasonable bar profile based on the Ruessink et
%  al, 2003, paper (Intersite comparisons of interannual nearshore bar
%  behavior, JGR).  The user must supply the background (mean) profile, the
%  a bar profile is superimposed.  The background profile must be
%  monotonic.
%  Data are output at all the supplied locations, x.  x must be regularly
%  spaced and IS ASSUMED TO START FROM [x,h] = [0,0]. Ensure that xb is
%  also corrected to shoreline-based coord system
%
% Inputs:
%   hShore, hSea    - landward and seaward limits of barred region (Ruessink)
%   x, hBackground  - user-supplied background (mean) beach profile
%   xb              - user-supplied location of a sand bar crest
% Output:
%   hBarred         - barred bathymetry (sum of background plus barred)

%  Holman, 10/11

% 'universal parameters from Ruessink et al, 03
del = 0.3;                  % min significant bar amplitude
La = 100;                   % for reconstructing bar length, L.
Lb = 0.27;
Sa = 0.53;                  % values for envelope amplitude, S
Sb = 0.57;
Sc = 0.09;

if (any(h0<0))           % ensure shore-based coord system
    error('Background profile depths must all be >= 0')
end
dx = mean(diff(x));              % for integration

% create the amplitude envelope, S.  Unlike Ruessink, we force S=0 at x=0
% and ramp up to S=del at hOff.  This also requires adjusting Smax to
% compensate for this ramp.  Seaward of hOff, we apply an artificial
% exponential taper, matching the value and slope at that location

SMax = 0.2*hSea;            % from Ruessink data
S0 = exp(-((1-((h0-hShore)/(hSea-hShore))).^Sa - Sb).^2/Sc); % env. shape
offInd = find(h0<=hSea, 1, 'last'); % find offshore limit of S-form
maxInd = S0(1:offInd)==max(S0(1:offInd));
SMaxCorr = SMax-x(maxInd)/x(offInd)*del;    % adjust so S=SMax with ramp
S = x/x(offInd)*del + SMax*exp(-((1-((h0-hShore)/(hSea-hShore))).^Sa - Sb).^2/Sc);

% add exponential taper to offshore end
dSdhOff = (S(offInd) - S(offInd-1))/(h0(offInd) - h0(offInd-1)); % slope
kExp = -dSdhOff/S(offInd); 
S(h0>hSea) = S(offInd)*exp(-kExp*(h0(h0>hSea)-h0(offInd)));

% now create the sand bar phase function
L = La*exp(Lb*h0);
kOfX = 2*pi./L;
theta = fliplr(cumsum(fliplr(kOfX)))*dx;   % integrate from offshore shoreward

% Now the full bar function.  First find psi from the bar position xb
hb = interp1(x,h0,xb);
thetab = interp1(h0, theta, hb);
psi = thetab;       % max at phase = 0;
hBar = -S .* cos(theta - psi);
h = h0 + hBar;



