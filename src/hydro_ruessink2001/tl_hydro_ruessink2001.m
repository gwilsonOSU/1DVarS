function [tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr]=tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tauw,tl_detady,tl_dgamma,bkgd)%,outvar)
%
% [tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr]=tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_detady,tl_dgamma,bkgd)
%
% TL-code for hydro_ruessink2001.m.  Background state 'bkgd' can be a struct taken
% directly from output of hydro_ruessink2001.m
%

[g,alpha,beta0,nu,rho,hmin,gammaType,betaType]=hydroParams();

% break out the bkgd vars
Ew   =bkgd.Ew;  % convert back from W/m2
Er   =bkgd.Er;
Db=bkgd.Db;
Dr=bkgd.Dr;
c    =bkgd.c    ;
cg   =bkgd.cg   ;
k    =bkgd.k    ;
h    =bkgd.h    ;
n    =bkgd.n    ;
theta=bkgd.theta;
omega=bkgd.omega;
x    =bkgd.x    ;
H    =bkgd.Hrms ;
dSxydx=bkgd.dSxydx;
Fy=bkgd.Fy;
v=bkgd.vbar;
H0=bkgd.H0;
theta0=bkgd.theta0;
gamma=bkgd.gamma;
dgamma=bkgd.dgamma;
Hm=bkgd.Hm;
Qb=bkgd.Qb;
ka_drag=bkgd.ka_drag;
detady=bkgd.detady;
beta=bkgd.beta;

h(h<hmin)=hmin;  % min depth constraint

nx=length(x);
dx=diff(x(1:2));

%-----------------------------
% begin tangent-linear code
%-----------------------------

if(strcmp(betaType,'const') | strcmp(betaType,'none'))
  tl_beta=zeros(nx,1);
end
if(strcmp(betaType,'none'))
  tl_Dr=zeros(nx,1);
end
if(~exist('dgamma'))
  tl_dgamma=zeros(nx,1);
end

tl_tauw=tl_tauw(:,2);

% dispersion
tl_k=-tl_h.*k.^2.*sech(k.*h).^2./(tanh(k.*h)+k.*h.*sech(k.*h).^2) ...
     + 2*omega/g./(tanh(k.*h) ...
     + k.*h.*sech(k.*h).^2).*tl_omega;
tl_c=-omega./k.^2.*tl_k + tl_omega./k;
% n= .5* + k.*h./sinh(2*k.*h);
tl_n = tl_k.*h./sinh(2*k.*h) + k.*tl_h./sinh(2*k.*h) ...
       - k.*h./sinh(2*k.*h).^2.*cosh(2*k.*h)*2.*(tl_k.*h+k.*tl_h);
tl_cg=tl_n.*c+n.*tl_c;
refconst=sin(theta0)/c(1);
tl_refconst=cos(theta0)/c(1)*tl_theta0-sin(theta0)/c(1)^2*tl_c(1);

% gamma can be either calculated based on deep water wave steepness (s0)
% following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
% or based on the empirical fit obtained for duck94 by Ruessink et
% al. (2003).
if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);  % = 2*pi*g/omega^2
  tl_L0 = -2*2*pi*g/omega^3*tl_omega;
  s0=H0/L0;
  tl_s0=tl_H0/L0-H0/L0^2*tl_L0;
  gamma=0.5+0.4*tanh(33*s0);
  tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0;
  gamma=ones(nx,1)*gamma;
  tl_gamma=ones(nx,1)*gamma;
elseif(gammaType==2003)
  gamma=0.76*k.*h+0.29;
  tl_gamma = 0.76*tl_k.*h + 0.76*k.*tl_h;
end
gamma=gamma+dgamma;
tl_gamma=tl_gamma+tl_dgamma;

% stepping, explicit scheme
tl_Ew=zeros(nx,1);
tl_Er=zeros(nx,1);
tl_Db=zeros(nx,1);
tl_Dr=zeros(nx,1);
tl_Ew(1)=rho*g/8*2*H0*tl_H0;
tl_Er(1)=0;
tl_H(1)=tl_H0;
tl_theta(1)=tl_theta0;
for i=2:nx

  % refraction
  % theta(i)=asin(c(i).*refconst);
  tl_theta(i)=1./sqrt(1-(c(i).*refconst).^2).*( refconst.*tl_c(i) + tl_refconst.*c(i) );

  % max wave height
  tharg=gamma(i-1)/0.88.*k(i-1).*h(i-1);
  tl_tharg=gamma(i-1)/0.88*( tl_k(i-1).*h(i-1) + k(i-1)*tl_h(i-1) ) ...
           + tl_gamma(i-1)/0.88.*k(i-1).*h(i-1);
  % Hm(i-1)=0.88./k(i-1).*tanh(tharg);
  tl_Hm(i-1)=0.88*( -1./k(i-1).^2.*tanh(tharg).*tl_k(i-1) ...
               + 1./k(i-1).*sech(tharg).^2.*tl_tharg );

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i-1)/Hm(i-1);
  tl_B=tl_H(i-1)/Hm(i-1)-H(i-1)/Hm(i-1)^2*tl_Hm(i-1);
  if(B<=.5)
    Qo=0;
    tl_Qo=0;
  else
    Qo=(2*B-1)^2;
    tl_Qo=2*(2*B-1).*2*tl_B;
  end
  if(B<=.2)
    tl_Qb(i-1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    tl_args=tl_Qo/B^2-2*(Qo-1)/B^3*tl_B;
    nums=Qo-exp(args);
    tl_nums=tl_Qo-exp(args)*tl_args;
    dens=B^2-exp(args);
    tl_dens=2*B*tl_B-exp(args)*tl_args;
    Qb(i-1)=Qo-B^2*nums/dens;
    tl_Qb(i-1)=tl_Qo-2*B*tl_B*nums/dens ...
        -B^2*tl_nums/dens + B^2*nums/dens^2*tl_dens;
  else
    tl_Qb(i-1)=0;
  end

  % term 1
  c1=alpha/4*rho*g*(omega/2/pi);  % alpha/4*rho*g/2/pi*omega
  tl_c1=alpha/4*rho*g/2/pi*tl_omega;
  % Db(i-1)=c1*Qb(i-1)*Hm(i-1)^2;
  tl_Db(i-1)=c1*tl_Qb(i-1)*Hm(i-1)^2 ...
      + 2*c1*Qb(i-1)*Hm(i-1)*tl_Hm(i-1) ...
      + tl_c1*Qb(i-1)*Hm(i-1)^2;

  % term 2
  nums1=cg(i-1)*Ew(i-1)*cos(theta(i-1));
  nums2=Db(i-1)*dx;
  denoms=cg(i)*cos(theta(i));
  tl_nums1=tl_cg(i-1)*Ew(i-1)*cos(theta(i-1)) ...
           + cg(i-1)*tl_Ew(i-1)*cos(theta(i-1)) ...
           - cg(i-1)*Ew(i-1)*sin(theta(i-1))*tl_theta(i-1);
  tl_nums2=tl_Db(i-1)*dx;
  tl_denoms=tl_cg(i)*cos(theta(i)) ...
            - cg(i)*sin(theta(i))*tl_theta(i);
  % Ew(i)=(nums1-nums2)/denoms;
  tl_Ew(i) = (tl_nums1-tl_nums2)/denoms ...
            - (nums1-nums2)/denoms^2*tl_denoms;

  % roller
  if(~strcmp(betaType,'none'))
    if(strcmp(betaType,'rafati21'))  % rafati et al. (2021) variable-beta
      if(k(i-1)*h(i-1)<0.45)
        tl_beta(i-1)=0;
      else
        if(beta(i-1)==0.1)  % limiter was hit
          tl_beta(i-1)=0;
        else
          tl_beta(i-1) = ...
              + 0.03*tl_k(i-1)*h(i-1)*(h(i-1)-H(i-1))/H(i-1) ...
              + 0.03*k(i-1)*tl_h(i-1)*(h(i-1)-H(i-1))/H(i-1) ...
              + 0.03*k(i-1)*h(i-1)*(tl_h(i-1)-tl_H(i-1))/H(i-1) ...
              - 0.03*k(i-1)*h(i-1)*(h(i-1)-H(i-1))/H(i-1)^2*tl_H(i-1);
        end
      end
    end
    nums=2*Er(i-1)*c(i-1)*cos(theta(i-1))+dx*(Db(i-1)-Dr(i-1));
    denoms=2*c(i)*cos(theta(i));
    tl_Dr(i-1) = ...
        + 2*g*tl_Er(i-1)*sin(beta(i-1))/c(i-1) ...
        + 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)*tl_beta(i-1) ...
        - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*tl_c(i-1);
    tl_nums = ...
        + 2*tl_Er(i-1)*c(i-1)*cos(theta(i-1)) ...
        + 2*Er(i-1)*tl_c(i-1)*cos(theta(i-1)) ...
        - 2*Er(i-1)*c(i-1)*sin(theta(i-1))*tl_theta(i-1) ...
        + dx*(tl_Db(i-1)-tl_Dr(i-1));
    tl_denoms = ...
        + 2*tl_c(i)*cos(theta(i)) ...
        - 2*c(i)*sin(theta(i))*tl_theta(i);
    if(Er(i)<0)
      tl_Er(i)=0;
    else
      tl_Er(i) = tl_nums/denoms - nums/denoms^2*tl_denoms;
    end
  end  % roller

  % H(i)=sqrt(8/rho/g*Ew(i));
  if(Ew(i)==0)
    tl_H(i)=0;
  else
    tl_H(i)=.5./sqrt(8/rho/g*Ew(i))*8/rho/g.*tl_Ew(i);
  end

end
tl_c=tl_c(:);
tl_H=tl_H(:);
tl_theta=tl_theta(:);

% radiation stress gradient
if(~strcmp(betaType,'none'))
  % dSxydx = -sin(theta)./c.*Dr/rho;
  tl_dSxydx = -cos(theta)./c.*Dr/rho.*tl_theta ...
      +sin(theta)./c.^2.*Dr/rho.*tl_c ...
      -sin(theta)./c.*tl_Dr/rho;
else
  % dSxydx = -sin(theta)./c.*Db/rho;
  tl_dSxydx = -cos(theta)./c.*Db/rho.*tl_theta ...
      +sin(theta)./c.^2.*Db/rho.*tl_c ...
      -sin(theta)./c.*tl_Db/rho;
end

% total force = radiation stress gradient + wind stress
% Fy=dSxydx+tauw/rho+g*h.*detady;
tl_Fy = tl_dSxydx + tl_tauw/rho + g*tl_h.*detady + g*h.*tl_detady;

% define mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=0; %[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000).  To get TL model, differentiate the eqn for v (i.e., the
% fsolve() line in waveModel.m) on both sides, then solve for
% tl_v
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
tl_Cd=0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*tl_h+tl_ka_drag./h);
urms=1.416*H*omega./(4*sinh(k.*h));
tl_urms=1.416*omega*( tl_H./(4*sinh(k.*h)) ...
                      -H./(4*sinh(k.*h).^2).*cosh(k.*h).*( tl_k.*h+k.*tl_h ) ) ...
        + 1.416*H./(4*sinh(k.*h))*tl_omega;
B=a^2+(v./urms).^2;
tl_N = tl_Fy./sqrt(B) + tl_urms.*(v.*Cd-v.^3.*Cd./(urms.^2.*B)) + tl_Cd.*v.*urms;
dens = -urms.*Cd - v.^2.*Cd./urms./B;
if(nu==0)
  tl_v=tl_N./dens;
else
  tl_v = inv(diag(dens)+A)*tl_N;  % tl_N = (dens + A) * tl_v
end

% % TEST: look at a specific variable
% tl_H=eval(['tl_' outvar]);
% tl_theta=zeros(nx,1);
% tl_v=zeros(nx,1);
% tl_k=zeros(nx,1);
% tl_Ew=zeros(nx,1);
% tl_Er=zeros(nx,1);
% tl_Dr=zeros(nx,1);
