function [H,theta,v,k,Ew,Er,Dr,bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tauw,detady,dgamma)
%
% [H,theta,v,k,Ew,Er,Dr,bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tauw,detady,dgamma)
%
% Wave energy balance equation solver, explicit spatial stepping scheme,
% followed by longshore current momentum balance solver with constant
% horizontal mixing.
%
% Breaking dissipation using TG1983
% Roller energy and roller dissipation following Reniers & Battjes (1996)
%
% INPUTS:
%
% x       : grid is +'ve onshore, and x(1) = offshore boundary
% h       : water depth, m
% H0      : rms wave height at offshore boundary, m
% theta0  : wave angle at offshore boundary, rads
% omega   : wave frequency, rad/m
% ka_drag : hydraulic roughness factor, m
% tauw    : vector (x,y) components of wind stress, N/m2.  Only tauw(:,2) (longshore) is actually used.
% detady  : longshore pressure gradient, m/m units
%
% OUTPUTS:
%
% H     : rms wave height, m
% theta : wave angle, rads
% v     : longshore current, m/s
% k     : wavenumber, rad/m
% wkspc : struct containing all NL variables, used for TL-AD models
%

[g,alpha,beta,nu,rho,hmin,gammaType]=hydroParams();

h(h<hmin)=hmin;  % min depth constraint

tauw2d=tauw;
tauw=tauw(:,2);  % for compatibility

% grid
nx=length(x);
dx=diff(x(1:2));

if(~exist('dgamma'))
  dgamma=zeros(nx,1);  % for compatibility
end

% dispersion
k=fsolve(@(k)omega^2-g*k.*tanh(k.*h),omega./sqrt(g*h),optimset('Display','off'));
c=max(0,real(omega./k));
n=.5*(1+2*k.*h./sinh(2*k.*h));
cg=n.*c;
refconst=sin(theta0)/c(1);

% gamma can be either calculated based on deep water wave steepness (s0)
% following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
% or based on the empirical fit obtained for duck94 by Ruessink et
% al. (2003).
if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);
  s0=H0/L0;
  gamma=0.5+0.4*tanh(33*s0);
  gamma=ones(nx,1)*gamma;
elseif(gammaType==2003)
  gamma=0.76*k.*h+0.29;
end
gamma=gamma+dgamma;

% refraction
theta=asin(c.*refconst);

% stepping, explicit scheme
Ew=zeros(nx,1);
Er=zeros(nx,1);
Db=zeros(nx,1);
Dr=zeros(nx,1);
Ew(1)=rho*g/8*H0^2;
Er(1)=0;
H(1)=H0;
theta(1)=theta0; %asin(c(1).*refconst);
for i=2:nx

  % max wave height
  tharg=gamma(i-1)/0.88.*k(i-1).*h(i-1);
  Hm(i-1)=0.88./k(i-1).*tanh(tharg);

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i-1)/Hm(i-1);
  if(B<=.5)
    Qo=0;
  else
    Qo=(2*B-1)^2;
  end
  if(B<=.2)
    Qb(i-1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    nums=Qo-exp(args);
    dens=B^2-exp(args);
    Qb(i-1)=Qo-B^2*nums/dens;
  else
    Qb(i-1)=1;
  end
  c1=alpha/4*rho*g*(omega/2/pi);
  Db(i-1)=c1*Qb(i-1)*Hm(i-1)^2;

  nums1=cg(i-1)*Ew(i-1)*cos(theta(i-1));
  nums2=Db(i-1)*dx;
  denoms=cg(i)*cos(theta(i));
  Ew(i)=(nums1-nums2)/denoms;

  if(beta>0)
    Dr(i-1)=2*g*Er(i-1)*sin(beta)/c(i-1);
    Er(i)=(2*Er(i-1)*c(i-1)*cos(theta(i-1))+dx*(Db(i-1)-Dr(i-1)))/(2*c(i)*cos(theta(i)));
    if(Er(i)<0)
      Er(i)=0;
    end
  end

  if(Ew(i)<.001)
    Ew(i)=.001;
  end
  H(i)=sqrt(8/rho/g*Ew(i));

end
c=c(:);
cg=cg(:);
k=k(:);
n=n(:);
theta=theta(:);
H=H(:);

% radiation stress gradient
if(beta>0)  % roller
  dSxydx = -sin(theta)./c.*Dr/rho;
else
  dSxydx=-sin(theta)./c.*Db/rho;
end
dSxydx(dSxydx==0)=1e-6;  % avoid singularity in TL model

% bottom stress model following Ruessink et al. (2001), Feddersen et al. (2000)

% total force = radiation stress gradient + wind stress + pressure gradient,
% m2/s2 units
Fy=dSxydx+tauw/rho+g*h.*detady;

% v1: analytical solution, no mixing
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
urms=1.416*H.*omega./(4*sinh(k.*h));
v2 = sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2;
v=sqrt(v2).*sign(-Fy);

% mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=0; %[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);

% v2: nonlinear solution with mixing
v0=v;
opt=optimset('Display','off');
v = fsolve(@(v)Fy + Cd.*urms.*v.*sqrt(a^2+(v./urms).^2) - A*v,v0,opt);
v=real(v);

% outputs struct
bkgd.x=x;
bkgd.h=h;
bkgd.H0     =H0     ;
bkgd.theta0 =theta0 ;
bkgd.omega  =omega  ;
bkgd.ka_drag=ka_drag;
bkgd.tauw   =tauw2d ;
bkgd.Ew   =Ew;
bkgd.Er   =Er;
bkgd.Db=Db;
bkgd.Dr=Dr;
bkgd.c=c;
bkgd.cg=cg;
bkgd.k=k;
bkgd.n=n;
bkgd.theta=theta;
bkgd.omega=omega;
bkgd.Hrms=H;
bkgd.gamma=gamma;
bkgd.dgamma=dgamma;
bkgd.Hm=Hm;
bkgd.Qb=Qb;
bkgd.dSxydx=dSxydx;
bkgd.Fy=Fy;
bkgd.vbar=real(v);
bkgd.detady=detady;

% optional, if user only wants the struct
if(nargout==1)
  H=bkgd;
end
