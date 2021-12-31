function [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,phs,bkgd)
%
% [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,phs,bkgd)
%
% TL-code for Uwave_ruessink2012.m
%
% note, 'bkgd' struct should be taken from last output of
% Uwave_ruessink2012.m, this contains the NL calculated variables

% Break out NL background vars. Note this hard-coded version runs
% significantly faster than if I use eval() to dynamically load all bkgd
% variables
Uw =bkgd.Uw ;
B  =bkgd.B  ;
psi=bkgd.psi;
phi=bkgd.phi;
b  =bkgd.b  ;
r  =bkgd.r  ;
f  =bkgd.f  ;
f1 =bkgd.f1 ;
f2 =bkgd.f2 ;
u  =bkgd.u  ;
Aw =bkgd.Aw ;
Sw =bkgd.Sw ;

%----------------------------
% begin TL code
%----------------------------

% convert to B and psi
% B = sqrt(Aw.^2 + Sw.^2);
tl_B = 1./sqrt(Aw.^2 + Sw.^2).*(Aw.*tl_Aw+Sw.*tl_Sw);
% psi = atan(Aw./Sw);
tl_psi = 1./(1+(Aw./Sw).^2).*( tl_Aw./Sw - Aw./Sw.^2.*tl_Sw );

% phi parameter, eqn (12)
tl_phi = -tl_psi;  % verified

% r parameter, 2-step inversion.  First solve for b as a function of B, then
% solve for r as a function of b.  TL code then follows based on derivation
% of partial derivatives from the 2 trancendental functions being solved,
% see ./Uwave_TLbB_deriv/*.jpg for derivations
tl_b = sqrt(2)/3*(1-b.^2).^(3/2).*tl_B;    % uses analytical db/db
tl_r = (1+f)./(1+r.^2./(f.*(1+f))).*tl_b;  % uses analytical dr/db

% Abreu et al. (2010) time series formula, eqn (4).  Note tl_phs==0
% f = sqrt(1-r.^2);
tl_f = .5./f.*( -2*r.*tl_r );
% f1 = sin(phs) + r.*sin(phi)./(1+f);
tl_f1 = + tl_r.*sin(phi)./(1+f) ...
        + r.*cos(phi)./(1+f).*tl_phi ...
        - r.*sin(phi)./(1+f).^2.*tl_f;
% f2 = 1 - r.*cos(phs+phi);
tl_f2 = -tl_r.*cos(phs+phi) ...
        + r.*sin(phs+phi).*tl_phi;
% u = Uw.*f.*f1./f2;
tl_u = tl_Uw.*f.*f1./f2 ...
       + Uw.*tl_f.*f1./f2 ...
       + Uw.*f.*tl_f1./f2 ...
       - Uw.*f.*f1./f2.^2.*tl_f2;
