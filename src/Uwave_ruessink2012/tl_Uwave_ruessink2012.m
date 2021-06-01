function [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(tl_Hmo,tl_k,tl_omega,tl_h,phs,bkgd)
%
% [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(tl_Hmo,tl_k,tl_omega,tl_h,phs,bkgd)
%
% TL-code for Uwave_ruessink2012.m
%
% note, 'bkgd' struct should be taken from last output of
% Uwave_ruessink2012.m, this contains the NL calculated variables

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
bra=2*(1-b.^2);  % not calculated in NL code

%----------------------------
% begin TL code
%----------------------------

% Ursell number, eqn (6)
tl_aw=tl_Hmo/2;
tl_Ur = 3/4*( (tl_aw.*k+aw.*tl_k)./(k.*h).^3 ...
              - 3*aw.*k./(k.*h).^4.*( tl_k.*h + k.*tl_h ) );

% wave velocity magnitude, linear theory, stated in text
tl_Hrms=tl_Hmo/1.4;
tl_Uw = omega/2.*( tl_Hrms./sinh(k.*h) ...
                   - Hrms./sinh(k.*h).^2.*cosh(k.*h).*(tl_k.*h+k.*tl_h) ) ...
        + tl_omega/2.*Hrms./sinh(k.*h);

% non-linearity param B, eqn (9).  Using parameter values quoted in text
tl_ee=-1./Ur./p4.*tl_Ur;
tl_dens=exp(ee).*tl_ee;
tl_B = -(p2-p1)./dens.^2.*tl_dens;

% non-linearity phase param psi, eqn (10)
tl_psi = pi/2*sech(p5./Ur.^p6).^2.*( -p5.*p6.*Ur.^(-p6-1).*tl_Ur );

% phi parameter, eqn (12)
tl_phi = -tl_psi;  % verified

% r parameter, 2-step inversion.  First solve for b as a function of B, then
% solve for r as a function of b.  TL code then follows based on derivation
% of partial derivatives from the 2 trancendental functions being solved,
% see ./Uwave_TLbB_deriv/*.jpg for derivations
tl_b = sqrt(bra)./(3-6*b.^2./bra).*tl_B;  % uses db/dB eqn.
tl_r = (1+f)./(1+r.^2./(f.*(1+f))).*tl_b;  % uses dr/db
tl_f = .5./f.*( -2*r.*tl_r );  % ok
% TODO: tl_b and tl_r both look to be a bit systematically biased in testing...

% Abreu et al. (2010) time series formula, eqn (4)
tl_f1 = tl_r.*sin(phi)./(1+f) ...
        + r.*cos(phi)./(1+f).*tl_phi ...
        - r.*sin(phi)./(1+f).^2.*tl_f; % ok
tl_f2 = -tl_r.*cos(phs+phi) ...
        + r.*sin(phs+phi).*tl_phi;
tl_u = tl_Uw.*f.*f1./f2 ...
       + Uw.*tl_f.*f1./f2 ...
       + Uw.*f.*tl_f1./f2 ...
       - Uw.*f.*f1./f2.^2.*tl_f2;  % ok
