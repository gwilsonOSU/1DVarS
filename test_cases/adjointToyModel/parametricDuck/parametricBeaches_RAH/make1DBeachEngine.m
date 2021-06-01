function h = make1DBeachEngine(x,xs,bs,xb,xd,hd,bd,hsea)
%   h = make1DBeachEngine(x,xShore,betaShore,xBar,xDeep,hDeep, ...
%                           betaDeep, hSea)
%
%  Compute a 1D parametric beach following Holman et al, 2014, C. Eng.
%  which follows Ruessink et al, 2003.  
%  All x locations are entered in user coordinates but are corrected to
%  local (shore = 0) coords for the calculations.  
%  Inputs:
%       x       - N by vector of x locations to find h
%       xShore  - shoreline (h = 0) location
%       betaShore - CLIMATOLOGICAL shoreline beach slope
%       xBar    - x location of bar crest
%       xDeep   - x location of some location seaward of the active bar zone
%       hDeep   - depth that this offshore location
%       betaDeep - bathymetric slope at this offshore location
%       hSea    - depth of seaward limit of bar envelope (see Ruessink)
%  Output:
%       h       - depth at x locations, nan if landward of shoreline

h = nan(size(x));                % default to nan
good = find(x>=xs);             % solve only for wet locations
h0 = makeBackgroundProfileComposite(x(good)-xs,xd-xs,hd,bs,bd);
hb = makeBarredBeach(0,hsea,x(good),h0,xb);
h(good) = hb;
