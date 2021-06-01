function [ad_d50,ad_d90,ad_h,ad_tanbeta,...
          ad_Hrms,ad_kabs,ad_omega,ad_theta,ad_ubar,...
          ad_Dr,ad_param] = ad_qtrans_soulsbyVanRijn(ad_qtot,bkgd) %,invar)
%
% AD-code for qtrans_soulsbyVanRijn.m
%

% break out NL background vars
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
psi=bkgd.psi;

physicalConstants;

%-------------------------------------
% begin ad code
%-------------------------------------

% init ad variables
nx=length(h);
ad_d50         =zeros(nx,1);
ad_d90         =zeros(nx,1);
ad_h           =zeros(nx,1);
ad_tanbeta     =zeros(nx,1);
ad_Hrms        =zeros(nx,1);
ad_kabs        =zeros(nx,1);
ad_theta       =zeros(nx,1);
ad_ubar        =zeros(nx,2);
ad_Dr          =zeros(nx,1);
ad_param.facua =0;
ad_param.alphab=0;
ad_q=zeros(nx,1);
ad_Dh=zeros(nx,1);
ad_C=zeros(nx,1);
ad_dCdx=zeros(nx,1);
ad_Ass=0;
ad_Asb=zeros(nx,1);
ad_Ufact=zeros(nx,1);
ad_slopeFact=zeros(nx,1);
ad_uAV=zeros(nx,1);
ad_uE2=zeros(nx,1);
ad_urms=zeros(nx,1);
ad_Cd=zeros(nx,1);
ad_ucr=zeros(nx,1);
ad_Dstar=0;
ad_VW=zeros(nx,1);
ad_Sk=zeros(nx,1);
ad_As=zeros(nx,1);
ad_dens=zeros(nx,1);
ad_nums=zeros(nx,1);
ad_earg=zeros(nx,1);
ad_psi=zeros(nx,1);
ad_Ur=zeros(nx,1);
ad_aw=zeros(nx,1);
ad_Hmo=zeros(nx,1);
ad_omega=0;

% % TEST: redefine input var
% eval(['ad_' invar '=ad_qtot;'])
% if(~strcmp(invar,'qtot'))
%   ad_qtot=zeros(nx,1);
% end

if(length(h)>2)
  %25b tl_qtot = tl_q ...
  %           + tl_Dh.*ddx_upwind(x,C) ...
  %           + Dh.*tl_dCdx;
  ad_q   =ad_q   +                  ad_qtot;
  ad_Dh  =ad_Dh  + ddx_upwind(x,C).*ad_qtot;
  ad_dCdx=ad_dCdx+ Dh             .*ad_qtot;
  ad_qtot=0;
  %25a tl_dCdx=tl_ddx_upwind(tl_C,x,C);
  ad_C=ad_C+ad_ddx_upwind(ad_dCdx,x,C);
  ad_dCdx=0;
else
  %25 tl_qtot=tl_q;
  ad_q=ad_q+ad_qtot;
  ad_qtot=0;
end

%24 tl_Dh(ind) = facDc*nuhfac*tl_h(ind).*(Dr(ind)/rho).^(1/3) ...
%         + (1/3)*facDc*nuhfac*h(ind).*(Dr(ind)/rho).^(1/3-1).*tl_Dr(ind)/rho;
ind=find(Dr>0);
ad_h(ind) =ad_h(ind) + facDc*nuhfac.*(Dr(ind)/rho).^(1/3)                   .*ad_Dh(ind);
ad_Dr(ind)=ad_Dr(ind)+ (1/3)*facDc*nuhfac*h(ind).*(Dr(ind)/rho).^(1/3-1)/rho.*ad_Dh(ind);
ad_Dh=0;

ind=find(Ufact>0);

%23 tl_C(ind) = (tl_Ass+tl_Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind) ...
%     + 2.4*(Ass+Asb(ind)).*Ufact(ind).^(1.4).*slopeFact(ind)./h(ind).*tl_Ufact(ind) ...
%     + (Ass+Asb(ind)).*Ufact(ind).^(2.4).*tl_slopeFact(ind)./h(ind) ...
%     - (Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind).^2.*tl_h(ind);
coef1=Ufact(ind).^(2.4).*slopeFact(ind)./h(ind);
coef2=2.4*(Ass+Asb(ind)).*Ufact(ind).^(1.4).*slopeFact(ind)./h(ind);
coef3=(Ass+Asb(ind)).*Ufact(ind).^(2.4)./h(ind);
coef4=(Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind).^2;
ad_Ass           =ad_Ass           + sum(coef1.*ad_C(ind));
ad_Asb(ind)      =ad_Asb(ind)      + coef1.*ad_C(ind);
ad_Ufact(ind)    =ad_Ufact(ind)    + coef2.*ad_C(ind);
ad_slopeFact(ind)=ad_slopeFact(ind)+ coef3.*ad_C(ind);
ad_h(ind)        =ad_h(ind)        - coef4.*ad_C(ind);
ad_C(ind)=0;

%22 tl_q(ind) = (tl_Ass+tl_Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind) ...
%     + (Ass+Asb(ind)).*tl_uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind) ...
%     + 2.4*(Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4-1).*slopeFact(ind).*tl_Ufact(ind) ...
%     + (Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*tl_slopeFact(ind);
coef1=uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind);
coef2=(Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind);
coef3=2.4*(Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4-1).*slopeFact(ind);
coef4=(Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4);
ad_Ass           =ad_Ass           + sum(coef1.*ad_q(ind));
ad_Asb(ind)      =ad_Asb(ind)      + coef1.*ad_q(ind);
ad_uAV(ind)      =ad_uAV(ind)      + coef2.*ad_q(ind);
ad_Ufact(ind)    =ad_Ufact(ind)    + coef3.*ad_q(ind);
ad_slopeFact(ind)=ad_slopeFact(ind)+ coef4.*ad_q(ind);
ad_q(ind)=0;

%21 tl_C=zeros(length(x),1);
%20 tl_q=zeros(length(x),1);
ad_C=zeros(length(x),1);
ad_q=zeros(length(x),1);

%19 tl_slopeFact = -tl_param.alphab*tanbeta ...
%     - param.alphab*tl_tanbeta;
ad_param.alphab=sum(ad_param.alphab- tanbeta.*ad_slopeFact);
ad_tanbeta     =ad_tanbeta - param.alphab.*ad_slopeFact;
ad_slopeFact=0;

%18 tl_Ufact = .5./sqrt(uE2 + .018*urms.^2./Cd).*( ...
%     tl_uE2 + 2*.018*urms./Cd.*tl_urms - .018*urms.^2./Cd.^2.*tl_Cd ) ...
%     - tl_ucr;
coef=.5./sqrt(uE2 + .018*urms.^2./Cd);
ad_uE2 =ad_uE2 + coef                     .*ad_Ufact;
ad_urms=ad_urms+ 2*.018*coef.*urms./Cd    .*ad_Ufact;
ad_Cd  =ad_Cd  - .018*coef.*urms.^2./Cd.^2.*ad_Ufact;
ad_ucr=ad_ucr - ad_Ufact;
ad_Ufact=0;

%17 tl_Ass = .012*tl_d50*Dstar^(-.6)/((s-1)*g*d50)^1.2 ...
%          - .6*.012*d50*Dstar^(-.6-1)/((s-1)*g*d50)^1.2.*tl_Dstar ...
%          - 1.2*.012*d50*Dstar^(-.6)*((s-1)*g*d50)^(-2.2)*(s-1)*g*tl_d50;
coef1=.012*Dstar^(-.6)/((s-1)*g*d50)^1.2;
coef2=1.2*.012*d50*Dstar^(-.6)*((s-1)*g*d50)^(-2.2)*(s-1)*g;
coef3=.6*.012*d50*Dstar^(-.6-1)/((s-1)*g*d50)^1.2;
ad_d50  =ad_d50  + (coef1-coef2)*ad_Ass;
ad_Dstar=ad_Dstar- coef3        *ad_Ass;
ad_Ass=0;

%16 tl_Asb = .005*tl_h.*(d50./h).^1.2/((s-1)*g*d50)^1.2 ...
%          + 1.2*.005*h.*(d50./h).^(1.2-1)/((s-1)*g*d50)^1.2.*( ...
%              tl_d50./h - d50./h.^2.*tl_h ) ...
%          - 1.2*.005*h.*(d50./h).^1.2.*((s-1)*g*d50)^(-2.2).*(s-1)*g*tl_d50;
coef1=.005.*(d50./h).^1.2/((s-1)*g*d50)^1.2;
coef2=1.2*.005*h.*(d50./h).^(1.2-1)/((s-1)*g*d50)^1.2;
coef3=1.2*.005*h.*(d50./h).^1.2.*((s-1)*g*d50)^(-2.2).*(s-1)*g;
ad_h = ad_h + (coef1-coef2.*d50./h.^2).*ad_Asb;
ad_d50 = ad_d50 + sum((coef2./h - coef3).*ad_Asb);
ad_Asb=0;

%15 tl_Dstar=(g*(s-1)/nu^2)^(1/3)*tl_d50;
ad_d50=ad_d50+(g*(s-1)/nu^2)^(1/3)*ad_Dstar;
ad_Dstar=0;

%14 tl_uAV = tl_VW.*cos(theta) ...
%          - VW.*sin(theta).*tl_theta ...
%          + tl_ubar(:,1);
ad_VW       =ad_VW       + cos(theta)    .*ad_uAV;
ad_theta    =ad_theta    - VW.*sin(theta).*ad_uAV;
ad_ubar(:,1)=ad_ubar(:,1)+                 ad_uAV;
ad_uAV=0;

%13 tl_VW = tl_param.facua*urms.*(Sk-As) ...
%         + param.facua*tl_urms.*(Sk-As) ...
%         + param.facua*urms.*(tl_Sk-tl_As);
ad_param.facua=ad_param.facua+sum(urms.*(Sk-As)    .*ad_VW);
ad_urms       =ad_urms       + param.facua.*(Sk-As).*ad_VW;
ad_Sk         =ad_Sk         + param.facua*urms    .*ad_VW;
ad_As         =ad_As         - param.facua*urms    .*ad_VW;
ad_VW=0;

%12 tl_As = -.79./dens.^2.*sin(psi).*tl_dens ...
%         +.79./dens.*cos(psi).*tl_psi;
ad_dens=ad_dens-.79./dens.^2.*sin(psi).*ad_As;
ad_psi =ad_psi +.79./dens.*cos(psi)   .*ad_As;
ad_As=0;

%11 tl_Sk = -.79./dens.^2.*cos(psi).*tl_dens ...
%         - .79./dens.*sin(psi).*tl_psi;
ad_dens=ad_dens- .79./dens.^2.*cos(psi).*ad_Sk;
ad_psi =ad_psi - .79./dens.*sin(psi)   .*ad_Sk;
ad_Sk=0;

%10 tl_dens = exp(earg).*tl_earg;
ad_earg = ad_earg + exp(earg).*ad_dens;
ad_dens=0;

%9 tl_earg = 1./Ur/.35.*tl_Ur;
ad_Ur = ad_Ur + 1./Ur/.35.*ad_earg;
ad_earg=0;

%8 tl_psi = -pi/2*sech(.64./Ur.^(.60)).^2.*.6*.64.*Ur.^(-.6-1).*tl_Ur;
ad_Ur = ad_Ur -pi/2*sech(.64./Ur.^(.60)).^2.*.6*.64.*Ur.^(-.6-1).*ad_psi;
ad_psi=0;

%7 tl_Ur = 3/4*tl_aw.*kabs./(kabs.*h).^3 ...
%         + 3/4*aw.*tl_kabs./(kabs.*h).^3 ...
%         - 3*3/4*aw.*kabs./(kabs.*h).^4.*( ...
%             tl_kabs.*h + kabs.*tl_h );
coef1=3/4.*kabs./(kabs.*h).^3;
coef2=3/4*aw./(kabs.*h).^3;
coef3=3*3/4*aw.*kabs./(kabs.*h).^4;
ad_aw  =ad_aw  + coef1      .*ad_Ur;
ad_kabs=ad_kabs+ coef2      .*ad_Ur;
ad_kabs=ad_kabs- coef3.*h   .*ad_Ur;
ad_h   =ad_h   - coef3.*kabs.*ad_Ur;
ad_Ur=0;

%6 tl_aw=tl_Hmo/2;
ad_Hmo = ad_Hmo + ad_aw/2;
ad_aw=0;

%5 tl_Hmo = 1.4*tl_Hrms;
ad_Hrms = ad_Hrms + 1.4*ad_Hmo;
ad_Hmo=0;

%4 tl_urms = tl_Hrms/2*omega./sinh(kabs.*h) ...
%           - Hrms/2*omega./sinh(kabs.*h).^2.*cosh(kabs.*h).*( ...
%               tl_kabs.*h + kabs.*tl_h ) ...
%          + Hrms/2./sinh(kabs.*h)*tl_omega;
coef=Hrms/2*omega./sinh(kabs.*h).^2.*cosh(kabs.*h);
ad_Hrms=ad_Hrms+ 1/2*omega./sinh(kabs.*h).*ad_urms;
ad_kabs=ad_kabs- coef.*h                 .*ad_urms;
ad_h   =ad_h   - coef.*kabs              .*ad_urms;
ad_omega=ad_omega+ sum(Hrms/2./sinh(kabs.*h).*ad_urms);
ad_urms=0;

%3 tl_ucr = .1*0.19*d50^(.1-1)*log10(4*h/d90).*tl_d50 ...
%          + 0.19*d50^.1*d90./(4*h*log(10)).*( ...
%              + 4*tl_h/d90 ...
%              - 4*h/d90^2*tl_d90 );
coef1=.1*0.19*d50^(.1-1)*log10(4*h./d90);
coef2=0.19*d50^.1.*d90./(4*h*log(10));
ad_d50=ad_d50+ coef1          .*ad_ucr;
ad_d90=ad_d90- 4*coef2.*h./d90.^2.*ad_ucr;
ad_h  =ad_h  + 4*coef2./d90    .*ad_ucr;
ad_ucr=0;

%2 tl_uE2 = 2*ubar(:,1).*tl_ubar(:,1) ...
%          + 2*ubar(:,2).*tl_ubar(:,2);
ad_ubar(:,1)=ad_ubar(:,1)+ 2*ubar(:,1).*ad_uE2;
ad_ubar(:,2)=ad_ubar(:,2)+ 2*ubar(:,2).*ad_uE2;
ad_uE2=0;

%1 tl_Cd = -2*.4^2./(log(h/z0)-1).^3.*tl_h./h;
ad_h =ad_h -2*.4^2./(log(h/z0)-1).^3.*ad_Cd./h;
ad_Cd=0;
