function h = findDepthOnTransect(t, hSea)
%   trans = findDepthOnTransect(trans, hSea);
%
% find a parametric beach profile for a transect in the structure trans
% then find the depth at distance t.d

xt = [0: max([t.bar.d t.d])+1];    % use min amount of profile to span bar
ht = makeBackgroundProfileComposite(xt,t.deep.d,t.deep.h,t.beta0,t.deep.beta);
hb = makeBarredBeach(0,hSea,xt,ht,t.bar.d);
h = interp1(xt,hb, t.d);
    