function [ad_Hmo,ad_k,ad_omega,ad_h]=ad_Uwave_ruessink2012(ad_u,ad_r,ad_phi,phs,bkgd)
%
% AD-code for tl_Uwave_ruessink2012.m
%

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
bra=2*(1-b.^2);  % not calculated in NL code

%----------------------------
% begin AD code
%----------------------------

% initialize ad variables
nt=length(phs);
ad_f1=zeros(1,nt);
ad_f2=zeros(1,nt);
ad_B=0;
ad_Ur=0;
ad_dens=0;
ad_ee=0;
ad_Hrms=0;
ad_k=0;
ad_h=0;
ad_Hmo=0;
ad_aw=0;
ad_psi=0;
ad_f=0;
ad_Uw=0;
ad_b=0;
ad_omega=0;

%b07 Abreu et al. (2010) time series formula, eqn (4)

%3 tl_u = tl_Uw.*f.*f1./f2 ...
%        + Uw.*tl_f.*f1./f2 ...
%        + Uw.*f.*tl_f1./f2 ...
%        - Uw.*f.*f1./f2.^2.*tl_f2;  % ok
%
% % NOTE: this is an example of a scalar tl variable 'tl_Uw' multiplying a
% % vector 'f'.  To treat this adjoint correctly, need to sum over the vector.
% % To see this, consider the non-vectorized version of the code below.  Note
% % this pattern appears several times in this code.
%
% for i=1:nt
%   % tl_u(i)=tl_Uw*f(i)*f1(i)/f2(i);
%   ad_Uw=ad_Uw+*f(i)*f1(i)/f2(i)*ad_u(i);
%   ad_u(i)=0;
% end
ad_Uw=ad_Uw+ sum(f.*f1./f2       .*ad_u);
ad_f =ad_f + sum(Uw.*f1./f2      .*ad_u);
ad_f1=ad_f1+ Uw.*f./f2       .*ad_u;
ad_f2=ad_f2- Uw.*f.*f1./f2.^2.*ad_u;
ad_u=0*ad_u;
%2 tl_f2 = -tl_r.*cos(phs+phi) ...
%         + r.*sin(phs+phi).*tl_phi;
ad_r  =ad_r  - sum(cos(phs+phi)   .*ad_f2);
ad_phi=ad_phi+ sum(r.*sin(phs+phi).*ad_f2);
ad_f2=0*ad_f2;
%1 tl_f1 = tl_r.*sin(phi)./(1+f) ...
%         + r.*cos(phi)./(1+f).*tl_phi ...
%         - r.*sin(phi)./(1+f).^2.*tl_f; % ok
ad_r  =ad_r  + sum(sin(phi)./(1+f)      .*ad_f1);
ad_phi=ad_phi+ sum(r.*cos(phi)./(1+f)   .*ad_f1);
ad_f  =ad_f  - sum(r.*sin(phi)./(1+f).^2.*ad_f1);
ad_f1=0*ad_f1;

%b06 r parameter, 2-step inversion.  First solve for b as a function of B, then
% solve for r as a function of b.  TL code then follows based on derivation
% of partial derivatives from the 2 trancendental functions being solved,
% see ./Uwave_TLbB_deriv/*.jpg for derivations

%3 tl_f = .5./f.*( -2*r.*tl_r );  % ok
ad_r=ad_r- 1./f.*r.*ad_f;
ad_f=0;
%2 tl_r = (1+f)./(1+r.^2./(f.*(1+f))).*tl_b;  % uses dr/db
ad_b=ad_b+(1+f)./(1+r.^2./(f.*(1+f))).*ad_r;
ad_r=0;
%1 tl_b = sqrt(bra)./(3-6*b.^2./bra).*tl_B;  % uses db/dB eqn.
ad_B=ad_B+ sqrt(bra)./(3-6*b.^2./bra).*ad_b;
ad_b=0;

%b05 phi parameter, eqn (12)

%1 tl_phi = -tl_psi;  % verified
ad_psi=ad_psi-ad_phi;
ad_phi=0;

%b04 non-linearity phase param psi, eqn (10)

%1 tl_psi = pi/2*sech(p5./Ur.^p6).^2.*( -p5.*p6.*Ur.^(-p6-1).*tl_Ur );
ad_Ur=ad_Ur- pi/2*sech(p5./Ur.^p6).^2.*p5.*p6.*Ur.^(-p6-1).*ad_psi;
ad_psi=0;

%b03 non-linearity param B, eqn (9).  Using parameter values quoted in text

%3 tl_B = -(p2-p1)./dens.^2.*tl_dens;
ad_dens=ad_dens -(p2-p1)./dens.^2.*ad_B;
ad_B=0;
%2 tl_dens=exp(ee).*tl_ee;
ad_ee=ad_ee+exp(ee).*ad_dens;
ad_dens=0;
%1 tl_ee=-1./Ur./p4.*tl_Ur;
ad_Ur=ad_Ur-1./Ur./p4.*ad_ee;
ad_ee=0;

%b02 wave velocity magnitude, linear theory, stated in text

%2 tl_Uw = omega/2.*( tl_Hrms./sinh(k.*h) ...
%                    - Hrms./sinh(k.*h).^2.*cosh(k.*h).*(tl_k.*h+k.*tl_h) ) ...
%         + 1/2.*Hrms./sinh(k.*h)*tl_omega;
ad_Hrms=ad_Hrms+ omega/2./sinh(k.*h)                        .*ad_Uw;
ad_k   =ad_k   - omega/2.*Hrms./sinh(k.*h).^2.*cosh(k.*h).*h.*ad_Uw;
ad_h   =ad_h   - omega/2.*Hrms./sinh(k.*h).^2.*cosh(k.*h).*k.*ad_Uw;
ad_omega=ad_omega+ sum(1/2.*Hrms./sinh(k.*h)*ad_Uw);
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
