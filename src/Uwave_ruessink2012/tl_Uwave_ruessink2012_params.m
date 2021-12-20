function [tl_Aw,tl_Sw,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_k,tl_omega,tl_h,bkgd)

% Break out NL background vars. Note this hard-coded version runs
% significantly faster than if I use eval() to dynamically load all bkgd
% variables
aw   =bkgd.aw   ;
Ur   =bkgd.Ur   ;
Hrms =bkgd.Hrms ;
Uw   =bkgd.Uw   ;
p1   =bkgd.p1   ;
p2   =bkgd.p2   ;
p3   =bkgd.p3   ;
p4   =bkgd.p4   ;
p5   =bkgd.p5   ;
p6   =bkgd.p6   ;
ee   =bkgd.ee   ;
dens =bkgd.dens ;
B0   =bkgd.B0   ;
psi0 =bkgd.psi0 ;
Hmo  =bkgd.Hmo  ;
k    =bkgd.k    ;
omega=bkgd.omega;
h    =bkgd.h    ;
Aw   =bkgd.Aw   ;
Sw   =bkgd.Sw   ;

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
tl_ee=-1./Ur./p4.*tl_Ur/log(10);
tl_dens=exp(ee).*tl_ee;
tl_B0 = -(p2-p1)./dens.^2.*tl_dens;

% non-linearity phase param psi, eqn (10)
tl_psi0 = pi/2*sech(p5./Ur.^p6).^2.*( -p5.*p6.*Ur.^(-p6-1).*tl_Ur );

% convert to Aw and Sw, including user provided corrections
% Aw = B0.*sin(psi0);
tl_Aw = + sin(psi0).*tl_B0 ...
        + B0.*cos(psi0).*tl_psi0;
% Sw = B0.*cos(psi0);
tl_Sw = + cos(psi0).*tl_B0 ...
        - B0.*sin(psi0).*tl_psi0;
