function [ad_Aw,ad_Sw,ad_Uw]=ad_Uwave_ruessink2012(ad_u,ad_r,ad_phi,phs,bkgd)
%
% AD-code for tl_Uwave_ruessink2012.m
%

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
ad_b=0;
ad_omega=0;
ad_Uw=0;
ad_Aw=0;
ad_Sw=0;

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
%0 tl_f = .5./f.*( -2*r.*tl_r );  % ok
ad_r=ad_r- 1./f.*r.*ad_f;
ad_f=0;

%b06 r parameter, 2-step inversion.  First solve for b as a function of B, then
% solve for r as a function of b.  TL code then follows based on derivation
% of partial derivatives from the 2 trancendental functions being solved,
% see ./Uwave_TLbB_deriv/*.jpg for derivations
%2 tl_r = (1+f)./(1+r.^2./(f.*(1+f))).*tl_b;  % uses dr/db
ad_b=ad_b+(1+f)./(1+r.^2./(f.*(1+f))).*ad_r;
ad_r=0;
% tl_b = sqrt(2)/3*(1-b.^2).^(3/2).*tl_B;    % uses analytical db/db
ad_B=ad_B+ sqrt(2)/3*(1-b.^2).^(3/2).*ad_b;
ad_b=0;

%b05 phi parameter, eqn (12)
%1 tl_phi = -tl_psi;  % verified
ad_psi=ad_psi-ad_phi;
ad_phi=0;

% convert to B and psi
%2 tl_psi = 1./(1+(Aw./Sw).^2).*( tl_Aw./Sw - Aw./Sw.^2.*tl_Sw );
ad_Aw=ad_Aw+ 1./(1+(Aw./Sw).^2)./Sw       .*ad_psi;
ad_Sw=ad_Sw- 1./(1+(Aw./Sw).^2).*Aw./Sw.^2.*ad_psi;
ad_psi=0;
%1 tl_B = 1./sqrt(Aw.^2 + Sw.^2).*(Aw.*tl_Aw+Sw.*tl_Sw);
ad_Aw=ad_Aw+ 1./sqrt(Aw.^2 + Sw.^2).*Aw.*ad_B;
ad_Sw=ad_Sw+ 1./sqrt(Aw.^2 + Sw.^2).*Sw.*ad_B;
ad_B=0;
