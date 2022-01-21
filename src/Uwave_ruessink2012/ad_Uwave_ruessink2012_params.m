function [ad_Hmo,ad_k,ad_omega,ad_h]=ad_Uwave_ruessink2012_params(ad_Aw,ad_Sw,ad_Uw,bkgd)

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

%----------------------------
% begin AD code
%----------------------------

% initialize ad variables
ad_B0=0;
ad_Ur=0;
ad_dens=0;
ad_ee=0;
ad_Hrms=0;
ad_k=0;
ad_h=0;
ad_Hmo=0;
ad_aw=0;
ad_psi0=0;
ad_omega=0;
ad_dSw=0;

% apply masking to points with zero wave height, avoids NaN
imask=find(Hmo==0);
% tl_Aw(imask)=0;
ad_Aw(imask)=0;
% tl_Sw(imask)=0;
ad_Sw(imask)=0;
% tl_Uw(imask)=0;
ad_Uw(imask)=0;

Sw(imask)=1;
Aw(imask)=1;
Uw(imask)=1;
Ur(imask)=1;
dens(imask)=1;
ee(imask)=1;

% convert to Aw and Sw, including user provided corrections
%2 tl_Sw = + cos(psi0).*tl_B0 ...
%         - B0.*sin(psi0).*tl_psi0;
ad_B0  =ad_B0  + cos(psi0)    .*ad_Sw;
ad_psi0=ad_psi0- B0.*sin(psi0).*ad_Sw;
ad_Sw=0;
%1 tl_Aw = + sin(psi0).*tl_B0 ...
%         + B0.*cos(psi0).*tl_psi0;
ad_B0  =ad_B0  + sin(psi0)    .*ad_Aw;
ad_psi0=ad_psi0+ B0.*cos(psi0).*ad_Aw;
ad_Aw=0;

%b04 non-linearity phase param psi, eqn (10)
%1 tl_psi0 = pi/2*sech(p5./Ur.^p6).^2.*( -p5.*p6.*Ur.^(-p6-1).*tl_Ur );
ad_Ur=ad_Ur- pi/2*sech(p5./Ur.^p6).^2.*p5.*p6.*Ur.^(-p6-1).*ad_psi0;
ad_psi0=0;

%b03 non-linearity param B, eqn (9).  Using parameter values quoted in text
%3 tl_B0 = -(p2-p1)./dens.^2.*tl_dens;
ad_dens=ad_dens -(p2-p1)./dens.^2.*ad_B0;
ad_B0=0;
%2 tl_dens=exp(ee).*tl_ee;
ad_ee=ad_ee+exp(ee).*ad_dens;
ad_dens=0;
%1 tl_ee=-1./Ur./p4.*tl_Ur/log(10);
ad_Ur=ad_Ur-1./Ur./p4/log(10).*ad_ee;
ad_ee=0;

%b02 wave velocity magnitude, linear theory, stated in text
%2 tl_Uw = omega/2.*( tl_Hrms./sinh(k.*h) ...
%                    - Hrms./sinh(k.*h).^2.*cosh(k.*h).*(tl_k.*h+k.*tl_h) ) ...
%         + 1/2.*Hrms./sinh(k.*h)*tl_omega;
ad_Hrms=ad_Hrms+ omega/2./sinh(k.*h)                        .*ad_Uw;
ad_k   =ad_k   - omega/2.*Hrms./sinh(k.*h).^2.*cosh(k.*h).*h.*ad_Uw;
ad_h   =ad_h   - omega/2.*Hrms./sinh(k.*h).^2.*cosh(k.*h).*k.*ad_Uw;
ad_omega=ad_omega+ nansum(1/2.*Hrms./sinh(k.*h).*ad_Uw);
ad_Uw=0;
%1 tl_Hrms=tl_Hmo/1.4;
ad_Hmo=ad_Hmo+ad_Hrms/1.4;
ad_Hrms=0;

%b01 Ursell number, eqn (6)
%2 tl_Ur = 3/4*( (tl_aw.*k+aw.*tl_k)./(k.*h).^3 ...
%               - 3*aw.*k./(k.*h).^4.*( tl_k.*h + k.*tl_h ) );
ad_aw=ad_aw+ 3/4.*k./(k.*h).^3        .*ad_Ur;
ad_k =ad_k + 3/4*aw./(k.*h).^3        .*ad_Ur;
ad_k =ad_k - 3/4*3*aw.*k./(k.*h).^4.*h.*ad_Ur;
ad_h =ad_h - 3/4*3*aw.*k./(k.*h).^4.*k.*ad_Ur;
ad_Ur=0;
%1 tl_aw=tl_Hmo/2;
ad_Hmo=ad_Hmo+ad_aw/2;
ad_aw=0;

% make sure H=0 points are not NaN, they should be ignored and set all i/o to zeros
imask=find(Hmo==0);
ad_Hmo(imask)=0;
ad_k(imask)=0;
ad_h(imask)=0;
