function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauw2d,ad_detady,ad_dgamma,ad_beta0]=ad_hydro_ruessink2001(ad_Hrms,ad_theta,ad_vbar,ad_kabs,ad_Ew,ad_Er,ad_Dr,bkgd)%,invar)
%
% AD-code for hydro_ruessink2001.m.  Background state 'bkgd' can be a struct taken
% directly from output of hydro_ruessink2001.m
%

% break out the bkgd vars
Ew   =bkgd.Ew;  % convert back from W/m2
Er   =bkgd.Er;
Db=bkgd.Db;
Dr=bkgd.Dr;
c    =bkgd.c    ;
cg   =bkgd.cg   ;
kabs    =bkgd.kabs    ;
h    =bkgd.h    ;
n    =bkgd.n    ;
theta=bkgd.theta;
omega=bkgd.omega;
x    =bkgd.x    ;
Hrms    =bkgd.Hrms ;
dSxydx=bkgd.dSxydx;
Fy=bkgd.Fy;
vbar=bkgd.vbar;
H0=bkgd.H0;
theta0=bkgd.theta0;
gamma=bkgd.gamma;
dgamma=bkgd.dgamma;
Hm=bkgd.Hm;
Qb=bkgd.Qb;
ka_drag=bkgd.ka_drag;
detady=bkgd.detady;
beta=bkgd.beta;
g        =bkgd.g        ;
alpha    =bkgd.alpha    ;
beta0    =bkgd.beta0    ;
nu       =bkgd.nu       ;
rho      =bkgd.rho      ;
hmin     =bkgd.hmin     ;
gammaType=bkgd.gammaType;
betaType =bkgd.betaType ;

h(h<hmin)=hmin;  % min depth constraint

nx=length(x);
dx=diff(x(1:2));

%-----------------------------
% begin AD code
%-----------------------------

% init ad vars
ad_h=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_ka_drag=0;
ad_tau_wind=zeros(nx,1);
ad_Cd=zeros(nx,1);
ad_urms=zeros(nx,1);
ad_Fy=zeros(nx,1);
ad_N=zeros(nx,1);
ad_dSxydx=zeros(nx,1);
ad_c=zeros(nx,1);
ad_nums=0;
ad_dens=0;
ad_nums1=0;
ad_nums2=0;
ad_denoms=0;
ad_args=0;
ad_Db=zeros(nx,1);
ad_cg=zeros(nx,1);
ad_Qb=zeros(nx,1);
ad_Qo=0;
ad_Hm=zeros(nx,1);
ad_B=0;
ad_tharg=0;
ad_gamma=zeros(nx,1);
ad_dgamma=zeros(nx,1);
ad_refconst=0;
ad_s0=0;
ad_n=zeros(nx,1);
ad_detady=zeros(nx,1);
ad_L0=0;
ad_c1=0;
ad_omega=0;
ad_beta=zeros(nx,1);
ad_Hmo=zeros(nx,1);
ad_beta0=0;

% % TEST: look at a specific variable
% ad_theta=zeros(nx,1);
% ad_vbar=zeros(nx,1);
% ad_kabs=zeros(nx,1);
% ad_Ew=zeros(nx,1);
% ad_Er=zeros(nx,1);
% ad_Dr=zeros(nx,1);
% % tl_Hrms=eval(['tl_' outvar]);
% eval(['ad_' invar '=ad_Hrms;']);
% if(~strcmp(invar,'Hrms'))
%   ad_Hrms=zeros(nx,1);
% end

% tl_Dr(1)=tl_Dr(2);
ad_Dr(2)=ad_Dr(2)+ad_Dr(1);
ad_Dr(1)=0;

%b7 % bottom stress model following Ruessink et al. (2001), Feddersen et
% % al. (2000).  To get TL model, differentiate the eqn for vbar (i.e., the
% % fsolve() line in waveModel.m) on both sides, then solve for tl_vbar
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
urms=1.416*Hrms*omega./(4*sinh(kabs.*h));
B=a^2+(vbar./urms).^2;
dens = -urms.*Cd - vbar.^2.*Cd./urms./B;

% mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=0; %[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);
if(nu==0)
  %4   tl_vbar=tl_N./dens;
  ad_N=ad_N+ad_vbar./dens;
  ad_vbar=0;
else
  %4   tl_vbar = inv(diag(dens)+A)*tl_N;
  ad_N = ad_N + inv(diag(dens)+A)'*ad_vbar;
  ad_vbar=0;
end
%3 tl_N = tl_Fy./sqrt(B) + tl_urms.*(vbar.*Cd-vbar.^3.*Cd./(urms.^2.*B)) + tl_Cd.*vbar.*urms;
ad_Fy  =ad_Fy  + 1./sqrt(B)                    .*ad_N;
ad_urms=ad_urms+ (vbar.*Cd-vbar.^3.*Cd./(urms.^2.*B)).*ad_N;
ad_Cd  =ad_Cd  + vbar.*urms                       .*ad_N;
ad_N=0;
%2 tl_urms=1.416*omega*( tl_Hrms./(4*sinh(kabs.*h)) ...
%                       -Hrms./(4*sinh(kabs.*h).^2).*cosh(kabs.*h).*( tl_kabs.*h+kabs.*tl_h ) ) ...
%         + 1.416*Hrms./(4*sinh(kabs.*h))*tl_omega;
ad_Hrms=ad_Hrms+ 1.416*omega./(4*sinh(kabs.*h))                    .*ad_urms;
ad_kabs=ad_kabs- 1.416*omega*Hrms./(4*sinh(kabs.*h).^2).*cosh(kabs.*h).*h.*ad_urms;
ad_h=ad_h- 1.416*omega*Hrms./(4*sinh(kabs.*h).^2).*cosh(kabs.*h).*kabs.*ad_urms;
ad_omega = ad_omega + sum(1.416*Hrms./(4*sinh(kabs.*h)).*ad_urms);
ad_urms=0;
%1 tl_Cd=0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*tl_h+tl_ka_drag./h);
coef=0.015*(1/3)*(ka_drag./h).^(-2/3);
ad_h      =ad_h      - coef.*ka_drag./h.^2.*ad_Cd;
ad_ka_drag=ad_ka_drag+ sum(coef./h            .*ad_Cd);
ad_Cd=0;

%b6 % total force = radiation stress gradient + wind stress
%1 tl_Fy = tl_dSxydx + tl_tau_wind/rho + g*tl_h.*detady + g*h.*tl_detady;
ad_dSxydx=ad_dSxydx+ad_Fy;
ad_tau_wind  =ad_tau_wind  +ad_Fy/rho;
ad_h = ad_h + g*detady.*ad_Fy;
ad_detady = ad_detady + g*h.*ad_Fy;
ad_Fy=0;

%b5 % radiation stress gradient
if(~strcmp(betaType,'none'))
  %1   tl_dSxydx = -cos(theta)./c.*Dr/rho.*tl_theta ...
  %       +sin(theta)./c.^2.*Dr/rho.*tl_c ...
  %       -sin(theta)./c.*tl_Dr/rho;
  ad_theta=ad_theta-cos(theta)./c.*Dr/rho   .*ad_dSxydx;
  ad_c    =ad_c    +sin(theta)./c.^2.*Dr/rho.*ad_dSxydx;
  ad_Dr   =ad_Dr   -sin(theta)./c/rho       .*ad_dSxydx;
  ad_dSxydx=0;
else
  %1   tl_dSxydx = -cos(theta)./c.*Db/rho.*tl_theta ...
  %       +sin(theta)./c.^2.*Db/rho.*tl_c ...
  %       -sin(theta)./c.*tl_Db/rho;
  ad_theta=ad_theta-cos(theta)./c.*Db/rho   .*ad_dSxydx;
  ad_c    =ad_c    +sin(theta)./c.^2.*Db/rho.*ad_dSxydx;
  ad_Db   =ad_Db   -sin(theta)./c/rho       .*ad_dSxydx;
  ad_dSxydx=0;
end

refconst=sin(theta0)/c(1);

% %b3 stepping, explicit scheme
for i=nx:-1:2

  % b3g
  if(Ew(i)==0)
    %   tl_Hrms(i)=0;
    ad_Hrms(i)=0;
  else
    %   tl_Hrms(i)=.5./sqrt(8/rho/g*Ew(i))*8/rho/g.*tl_Ew(i);
    ad_Ew(i)=ad_Ew(i)+.5./sqrt(8/rho/g*Ew(i))*8/rho/g.*ad_Hrms(i);
    ad_Hrms(i)=0;
  end

  % roller
  if(~strcmp(betaType,'none'))
    nums=2*Er(i-1)*c(i-1)*cos(theta(i-1))+dx*(Db(i-1)-Dr(i-1));
    denoms=2*c(i)*cos(theta(i));
    ad_nums=0;  % init AD for this gridpoint
    ad_denoms=0;  % init AD for this gridpoint
    if(Er(i)<0)
      % r5a tl_Er(i)=0;
      ad_Er(i)=0;
    else
      % r5b tl_Er(i) = tl_nums/denoms - nums/denoms^2*tl_denoms;
      ad_nums  =ad_nums  + 1/denoms     *ad_Er(i);
      ad_denoms=ad_denoms- nums/denoms^2*ad_Er(i);
      ad_Er(i)=0;
    end
    %r4 tl_denoms = ...
    %     + 2*tl_c(i)*cos(theta(i)) ...
    %     - 2*c(i)*sin(theta(i))*tl_theta(i);
    ad_c(i)    =ad_c(i)    + 2*cos(theta(i))     *ad_denoms;
    ad_theta(i)=ad_theta(i)- 2*c(i)*sin(theta(i))*ad_denoms;
    ad_denoms=0;
    %r3 tl_nums = ...
    %     + 2*tl_Er(i-1)*c(i-1)*cos(theta(i-1)) ...
    %     + 2*Er(i-1)*tl_c(i-1)*cos(theta(i-1)) ...
    %     - 2*Er(i-1)*c(i-1)*sin(theta(i-1))*tl_theta(i-1) ...
    %     + dx*(tl_Db(i-1)-tl_Dr(i-1));
    ad_Er(i-1)   =ad_Er(i-1)   + 2*c(i-1)*cos(theta(i-1))        *ad_nums;
    ad_c(i-1)    =ad_c(i-1)    + 2*Er(i-1)*cos(theta(i-1))       *ad_nums;
    ad_theta(i-1)=ad_theta(i-1)- 2*Er(i-1)*c(i-1)*sin(theta(i-1))*ad_nums;
    ad_Db(i-1)   =ad_Db(i-1)   + dx                              *ad_nums;
    ad_Dr(i-1)   =ad_Dr(i-1)   - dx                              *ad_nums;
    ad_nums=0;
    %r2 tl_Dr(i-1) = ...
    %     + 2*g*tl_Er(i-1)*sin(beta(i-1))/c(i-1) ...
    %     + 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)*tl_beta(i-1) ...
    %     - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*tl_c(i-1);
    ad_Er(i-1)  =ad_Er(i-1)  + 2*g*sin(beta(i-1))/c(i-1)          *ad_Dr(i-1);
    ad_beta(i-1)=ad_beta(i-1)+ 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)  *ad_Dr(i-1);
    ad_c(i-1)   =ad_c(i-1)   - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*ad_Dr(i-1);
    ad_Dr(i-1)=0;
    if(strcmp(betaType,'rafati21'))  % rafati et al. (2021) variable-beta
      if(kabs(i-1)*h(i-1)<0.45)
        %r1a tl_beta(i-1)=0;
        ad_beta(i-1)=0;
      else
        if(beta(i-1)==0.1)  % limiter was hit
          % tl_beta(i-1)=0;
          ad_beta(i-1)=0;
        else
          %r1b tl_beta(i-1) = ...
          %     + 0.03*tl_kabs(i-1)*h(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1) ...
          %     + 0.03*kabs(i-1)*tl_h(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1) ...
          %     + 0.03*kabs(i-1)*h(i-1)*(tl_h(i-1)-tl_Hrms(i-1))/Hrms(i-1) ...
          %     - 0.03*kabs(i-1)*h(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1)^2*tl_Hrms(i-1);
          ad_kabs(i-1)=ad_kabs(i-1)+ 0.03*h(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1)         *ad_beta(i-1);
          ad_h(i-1)=ad_h(i-1)+ 0.03*kabs(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1)         *ad_beta(i-1);
          ad_h(i-1)=ad_h(i-1)+ 0.03*kabs(i-1)*h(i-1)/Hrms(i-1)                  *ad_beta(i-1);
          ad_Hrms(i-1)=ad_Hrms(i-1)- 0.03*kabs(i-1)*h(i-1)/Hrms(i-1)                  *ad_beta(i-1);
          ad_Hrms(i-1)=ad_Hrms(i-1)- 0.03*kabs(i-1)*h(i-1)*(h(i-1)-Hrms(i-1))/Hrms(i-1)^2*ad_beta(i-1);
          ad_beta(i-1)=0;
        end
      end
    end
  end  % roller

  %b3e % term 2
  nums1=cg(i-1)*Ew(i-1)*cos(theta(i-1));
  nums2=Db(i-1)*dx;
  denoms=cg(i)*cos(theta(i));
  %4 tl_Ew(i) = (tl_nums1-tl_nums2)/denoms ...
  %           - (nums1-nums2)/denoms^2*tl_denoms;
  ad_nums1 =ad_nums1 + 1/denoms              *ad_Ew(i);
  ad_nums2 =ad_nums2 - 1/denoms              *ad_Ew(i);
  ad_denoms=ad_denoms- (nums1-nums2)/denoms^2*ad_Ew(i);
  ad_Ew(i)=0;
  %3 tl_denoms=tl_cg(i)*cos(theta(i)) ...
  %           - cg(i)*sin(theta(i))*tl_theta(i);
  ad_cg(i)   =ad_cg(i)   + cos(theta(i))      *ad_denoms;
  ad_theta(i)=ad_theta(i)- cg(i)*sin(theta(i))*ad_denoms;
  ad_denoms=0;
  %2 tl_nums2=tl_Db(i-1)*dx;
  ad_Db(i-1)=ad_Db(i-1)+ad_nums2*dx;
  ad_nums2=0;
  %1 tl_nums1=tl_cg(i-1)*Ew(i-1)*cos(theta(i-1)) ...
  %          + cg(i-1)*tl_Ew(i-1)*cos(theta(i-1)) ...
  %          - cg(i-1)*Ew(i-1)*sin(theta(i-1))*tl_theta(i-1);
  ad_cg(i-1)   =ad_cg(i-1)   + Ew(i-1)*cos(theta(i-1))        *ad_nums1;
  ad_Ew(i-1)   =ad_Ew(i-1)   + cg(i-1)*cos(theta(i-1))        *ad_nums1;
  ad_theta(i-1)=ad_theta(i-1)- cg(i-1)*Ew(i-1)*sin(theta(i-1))*ad_nums1;
  ad_nums1=0;

  %b3d % term 1
  c1=alpha/4*rho*g*(omega/2/pi);
  %1 tl_Db(i-1)=c1*tl_Qb(i-1)*Hm(i-1)^2 ...
  %     + 2*c1*Qb(i-1)*Hm(i-1)*tl_Hm(i-1) ...
  %     + tl_c1*Qb(i-1)*Hm(i-1)^2;
  ad_Qb(i-1)=ad_Qb(i-1)+ c1*Hm(i-1)^2        *ad_Db(i-1);
  ad_Hm(i-1)=ad_Hm(i-1)+ 2*c1*Qb(i-1)*Hm(i-1)*ad_Db(i-1);
  ad_c1 = ad_c1 + ad_Db(i-1)*Qb(i-1)*Hm(i-1)^2;
  ad_Db(i-1)=0;
  % tl_c1=alpha/4*rho*g/2/pi*tl_omega;
  ad_omega = ad_omega+alpha/4*rho*g/2/pi*ad_c1;
  ad_c1=0;

  %b3c % fraction of breaking waves, non-implicit approximation from SWAN code
  B=Hrms(i-1)/Hm(i-1);
  if(B<=.5)
    Qo=0;
  else
    Qo=(2*B-1)^2;
  end
  if(B<=.2)
    %3a   tl_Qb(i-1)=0;
    ad_Qb(i-1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    nums=Qo-exp(args);
    dens=B^2-exp(args);
    Qb(i-1)=Qo-B^2*nums/dens;
    %3d   tl_Qb(i-1)=tl_Qo-2*B*tl_B*nums/dens ...
    %       -B^2*tl_nums/dens + B^2*nums/dens^2*tl_dens;
    ad_Qo  =ad_Qo  +                 ad_Qb(i-1);
    ad_B   =ad_B   - 2*B*nums/dens  *ad_Qb(i-1);
    ad_nums=ad_nums- B^2/dens       *ad_Qb(i-1);
    ad_dens=ad_dens+ B^2*nums/dens^2*ad_Qb(i-1);
    ad_Qb(i-1)=0;
    %3c   tl_dens=2*B*tl_B-exp(args)*tl_args;
    ad_B   =ad_B   + 2*B      *ad_dens;
    ad_args=ad_args- exp(args)*ad_dens;
    ad_dens=0;
    %3b   tl_nums=tl_Qo-exp(args)*tl_args;
    ad_Qo  =ad_Qo  +           ad_nums;
    ad_args=ad_args- exp(args)*ad_nums;
    ad_nums=0;
    %3a   tl_args=tl_Qo/B^2-2*(Qo-1)/B^3*tl_B;
    ad_Qo=ad_Qo+ 1/B^2       *ad_args;
    ad_B =ad_B - 2*(Qo-1)/B^3*ad_args;
    ad_args=0;
  else
    %3a   tl_Qb(i-1)=0;
    ad_Qb(i-1)=0;
  end
  if(B<=.5)
    %2a   tl_Qo=0;
    ad_Qo=0;
  else
    %2a   tl_Qo=2*(2*B-1).*2*tl_B;
    ad_B=ad_B+2*(2*B-1).*2*ad_Qo;
    ad_Qo=0;
  end
  %1 tl_B=tl_Hrms(i-1)/Hm(i-1)-Hrms(i-1)/Hm(i-1)^2*tl_Hm(i-1);
  ad_Hrms(i-1) =ad_Hrms(i-1) + 1/Hm(i-1)       *ad_B;
  ad_Hm(i-1)=ad_Hm(i-1)- Hrms(i-1)/Hm(i-1)^2*ad_B;
  ad_B=0;

  %b3b % max wave height
  tharg=gamma(i-1)/0.88.*kabs(i-1).*h(i-1);
  %2 tl_Hm(i-1)=0.88*( -1./kabs(i-1).^2.*tanh(tharg).*tl_kabs(i-1) ...
  %              + 1./kabs(i-1).*sech(tharg).^2.*tl_tharg );
  ad_kabs(i-1)=ad_kabs(i-1)- 0.88./kabs(i-1).^2.*tanh(tharg).*ad_Hm(i-1);
  ad_tharg =ad_tharg + 0.88./kabs(i-1).*sech(tharg).^2.*ad_Hm(i-1);
  ad_Hm(i-1)=0;
  %1 tl_tharg=gamma(i-1)/0.88*( tl_kabs(i-1).*h(i-1) + kabs(i-1)*tl_h(i-1) ) ...
  %          + tl_gamma(i-1)/0.88.*kabs(i-1).*h(i-1);
  ad_kabs(i-1)=ad_kabs(i-1)+ gamma(i-1)/0.88*h(i-1)     *ad_tharg;
  ad_h(i-1)=ad_h(i-1)+ gamma(i-1)/0.88*kabs(i-1)     *ad_tharg;
  ad_gamma(i-1) =ad_gamma(i-1) + 1/0.88.*kabs(i-1).*h(i-1)*ad_tharg;
  ad_tharg=0;

  %b3a % refraction
  %1 tl_theta(i)=1./sqrt(1-(c(i).*refconst).^2).*( refconst.*tl_c(i) + tl_refconst.*c(i) );
  coef=1./sqrt(1-(c(i).*refconst).^2);
  ad_c(i)    =ad_c(i)    + coef.*refconst.*ad_theta(i);
  ad_refconst=ad_refconst+ coef.*c(i)    .*ad_theta(i);
  ad_theta(i)=0;

end  % end of stepping scheme loop

%b3? init variables for stepping scheme
%8 tl_theta(1)=tl_theta0;
ad_theta0=ad_theta0+ad_theta(1);
ad_theta(1)=0;
%7 tl_Hrms(1)=tl_H0;
ad_H0=ad_H0+ad_Hrms(1);
ad_Hrms(1)=0;
%6 tl_Er(1)=0;
ad_Er(1)=0;
%5 tl_Ew(1)=rho*g/8*2*H0*tl_H0;
ad_H0=ad_H0+rho*g/8*2*H0*ad_Ew(1);
ad_Ew(1)=0;
%4 tl_Dr=zeros(nx,1);
ad_Dr=zeros(nx,1);
%3 tl_Db=zeros(nx,1);
ad_Db=zeros(nx,1);
%2 tl_Er=zeros(nx,1);
ad_Er=zeros(nx,1);
%1 tl_Ew=zeros(nx,1);
ad_Ew=zeros(nx,1);

%b2 gamma can be either calculated based on deep water wave steepness (s0)
% following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
% or based on the empirical fit obtained for duck94 by Ruessink et
% al. (2003).
% tl_gamma=tl_gamma+tl_dgamma;
ad_dgamma=ad_dgamma+ad_gamma;  % do not clear ad_gamma b/c it's an incremental-add
if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);
  s0=H0/L0;
  %2 tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0;
  ad_s0=ad_s0+0.4*sech(33*s0).^2.*33*sum(ad_gamma);
  ad_gamma=0;
  %1 tl_s0=tl_H0/L0-H0/L0^2*tl_L0;
  ad_H0=ad_H0+ad_s0/L0;
  ad_L0=ad_L0-H0/L0^2*ad_s0;
  ad_s0=0;
  %0 tl_L0 = -2*2*pi*g/omega^3*tl_omega;
  ad_omega = ad_omega-2*2*pi*g/omega^3*ad_L0;
  ad_L0=0;
elseif(gammaType==2003)
  % tl_gamma = 0.76*tl_kabs.*h + 0.76*kabs.*tl_h;
  ad_kabs = ad_kabs + 0.76*h.*ad_gamma;
  ad_h = ad_h + 0.76*kabs.*ad_gamma;
  ad_gamma=0;
end

%b1 % dispersion
%5 tl_refconst=cos(theta0)/c(1)*tl_theta0-sin(theta0)/c(1)^2*tl_c(1);
ad_theta0=ad_theta0+ cos(theta0)/c(1)  *ad_refconst;
ad_c(1)  =ad_c(1)  - sin(theta0)/c(1)^2*ad_refconst;
ad_refconst=0;
%4 tl_cg=tl_n.*c+n.*tl_c;
ad_n=ad_n+ c.*ad_cg;
ad_c=ad_c+ n.*ad_cg;
ad_cg=0;
%3 tl_n = tl_kabs.*h./sinh(2*kabs.*h) + kabs.*tl_h./sinh(2*kabs.*h) ...
%        - kabs.*h./sinh(2*kabs.*h).^2.*cosh(2*kabs.*h)*2.*(tl_kabs.*h+kabs.*tl_h);
coef=kabs.*h./sinh(2*kabs.*h).^2.*cosh(2*kabs.*h)*2;
ad_kabs=ad_kabs+ (h./sinh(2*kabs.*h) - coef.*h).*ad_n;
ad_h=ad_h+ (kabs./sinh(2*kabs.*h) - coef.*kabs).*ad_n;
ad_n=0;
%2 tl_c=-omega./kabs.^2.*tl_kabs + tl_omega./kabs;
ad_kabs=ad_kabs-omega./kabs.^2.*ad_c;
ad_omega = ad_omega + sum(ad_c./kabs);
ad_c=0;
%1 tl_kabs=-tl_h.*kabs.^2.*sech(kabs.*h).^2./(tanh(kabs.*h)+kabs.*h.*sech(kabs.*h).^2) ...
%     + 2*omega/g./(tanh(kabs.*h)+kabs.*h.*sech(kabs.*h).^2)*tl_omega;
ad_h=ad_h-ad_kabs.*kabs.^2.*sech(kabs.*h).^2./(tanh(kabs.*h)+kabs.*h.*sech(kabs.*h).^2);
ad_omega = ad_omega + sum(2*omega/g./(tanh(kabs.*h)+kabs.*h.*sech(kabs.*h).^2).*ad_kabs);
ad_kabs=0;

% tl_tau_wind=tl_tauw2d(:,2);  % for compatibility reasons
ad_tauw2d(:,2)=ad_tau_wind;
ad_tau_wind=zeros(nx,1);

if(strcmp(betaType,'const'))
  % tl_beta=tl_beta0*ones(nx,1);
  ad_beta0=sum(ad_beta);
  ad_beta=0;
elseif(strcmp(betaType,'none'))
  % tl_beta=0;
  ad_beta=0;
  ad_Dr=zeros(nx,1);
end
if(~exist('dgamma'))
  ad_dgamma=zeros(nx,1);
end
