function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauwin,ad_detady,ad_dgamma]=ad_hydro_ruessink2001(ad_H,ad_theta,ad_v,ad_k,ad_Ew,ad_Er,ad_Dr,bkgd)

[g,alpha,beta0,nu,rho,hmin,gammaType,betaType]=hydroParams();

% break out the bkgd vars
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} '=bkgd.' fld{i} ';']);
end
beta=bkgd.beta;
gamma=bkgd.gamma;
H=Hrms;
v=vbar;

h(h<hmin)=hmin;  % min depth constraint

nx=length(x);
dx=diff(x(1:2));

%---------------------------------------------
% init AD
%---------------------------------------------

ad_h=zeros(nx,1);
ad_htot=zeros(nx,1);
ad_eta=zeros(nx,1);
ad_Sxx=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_ka_drag=0;
ad_tauw=zeros(nx,1);
ad_tauwin=zeros(nx,2);
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
ad_s0=0;
ad_n=zeros(nx,1);
ad_detady=zeros(nx,1);
ad_L0=0;
ad_c1=0;
ad_omega=0;
ad_beta=zeros(nx,1);

%---------------------------------------------
% begin AD
%---------------------------------------------

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000).  To get TL model, differentiate the eqn for v (i.e., the
% fsolve() line in waveModel.m) on both sides, then solve for tl_v
B=a^2+(v./urms).^2;
dens = -urms.*Cd - v.^2.*Cd./urms./B;
if(nu==0)
  %67a1 tl_v=tl_N./dens;
  ad_N = ad_N + ad_v./dens;
  ad_v=0;
else
  %67b1 tl_v = inv(diag(dens)+A)*tl_N;
  ad_N = ad_N + inv(diag(dens)+A)'*ad_v;
  ad_v=0;
end
%66 tl_N = tl_Fy./sqrt(B) + tl_urms.*(v.*Cd-v.^3.*Cd./(urms.^2.*B)) + tl_Cd.*v.*urms;
ad_Fy  =ad_Fy  + 1./sqrt(B)                    .*ad_N;
ad_urms=ad_urms+ (v.*Cd-v.^3.*Cd./(urms.^2.*B)).*ad_N;
ad_Cd  =ad_Cd  + v.*urms                       .*ad_N;
ad_N=0;
%65 tl_urms=1.416*omega*( tl_H./(4*sinh(k.*htot)) ...
%                       -H./(4*sinh(k.*htot).^2).*cosh(k.*htot).*( tl_k.*htot+k.*tl_htot ) ) ...
%         + 1.416*H./(4*sinh(k.*htot))*tl_omega;
ad_H=ad_H+ 1.416*omega./(4*sinh(k.*htot))                 .*ad_urms;
ad_k=ad_k- 1.416*omega*H./(4*sinh(k.*htot).^2).*cosh(k.*htot).*htot.*ad_urms;
ad_htot=ad_htot- 1.416*omega*H./(4*sinh(k.*htot).^2).*cosh(k.*htot).*k.*ad_urms;
ad_omega = ad_omega + sum(1.416*H./(4*sinh(k.*htot)).*ad_urms);
ad_urms=0;
%64 tl_Cd=0.015*(1/3)*(ka_drag./htot).^(-2/3).*(-ka_drag./htot.^2.*tl_htot+tl_ka_drag./htot);
coef=0.015*(1/3)*(ka_drag./htot).^(-2/3);
ad_hhtot  =ad_htot   - coef.*ka_drag./htot.^2.*ad_Cd;
ad_ka_drag=ad_ka_drag+ sum(coef./htot        .*ad_Cd);
ad_Cd=0;

% total force = radiation stress gradient + wind stress + pressure gradient,
% m2/s2 units
%63 tl_Fy = ...
%     + tl_dSxydx ...
%     + tl_tauw/rho ...
%     + g*tl_htot.*detady ...
%     + g*htot.*tl_detady;
ad_dSxydx=ad_dSxydx+ad_Fy;
ad_tauw  =ad_tauw  +ad_Fy/rho;
ad_h = ad_h + g*detady.*ad_Fy;
ad_detady = ad_detady + g*h.*ad_Fy;
ad_Fy=0;

% radiation stress gradient
if(~strcmp(betaType,'none'))
  %62a1 tl_dSxydx = ...
  %     - cos(theta)./c.*Dr/rho.*tl_theta ...
  %     + sin(theta)./c.^2.*Dr/rho.*tl_c ...
  %     - sin(theta)./c/rho.*tl_Dr;
  ad_theta=ad_theta-cos(theta)./c.*Dr/rho   .*ad_dSxydx;
  ad_c    =ad_c    +sin(theta)./c.^2.*Dr/rho.*ad_dSxydx;
  ad_Dr   =ad_Dr   -sin(theta)./c/rho       .*ad_dSxydx;
  ad_dSxydx=0;
else
  %62b1 tl_dSxydx = ...
  %     - cos(theta)./c.*Db/rho.*tl_theta ...
  %     + sin(theta)./c.^2.*Db/rho.*tl_c ...
  %     - sin(theta)./c/rho.*tl_Db;
  ad_theta=ad_theta-cos(theta)./c.*Db/rho   .*ad_dSxydx;
  ad_c    =ad_c    +sin(theta)./c.^2.*Db/rho.*ad_dSxydx;
  ad_Db   =ad_Db   -sin(theta)./c/rho       .*ad_dSxydx;
  ad_dSxydx=0;
end

% explicit forward stepping for wave propagation
for i=nx:-1:2
  khtot=k(i-1)*htot(i-1);

  % init AD constants for this loop iteration
  ad_term1=0;
  ad_term2=0;
  ad_khtot=0;
  ad_refconst=0;

  % update wave angle
  %61 tl_theta(i) = 1/sqrt(1-(c(i)*refconst)^2)*( refconst*tl_c(i) + tl_refconst*c(i) );
  ad_c(i)    =ad_c(i)    + 1/sqrt(1-(c(i)*refconst)^2)*refconst*ad_theta(i);
  ad_refconst=ad_refconst+ 1/sqrt(1-(c(i)*refconst)^2)*c(i)    *ad_theta(i);
  ad_theta(i)=0;

  % update wave induced setup
  term1=(cos(theta(i))^2+1)*cg(i)/c(i)-.5;
  term2=cos(theta(i))^2;
  %60 tl_htot(i) = tl_h(i) + tl_eta(i);
  ad_h(i)  =ad_h(i)  + ad_htot(i);
  ad_eta(i)=ad_eta(i)+ ad_htot(i);
  ad_htot(i)=0;
  %59 tl_eta(i) = ...
  %     + tl_eta(i-1) ...
  %     + (tl_Sxx(i-1)-tl_Sxx(i))/g/h(i-1)/rho ...
  %     - (Sxx(i-1)-Sxx(i))/g/h(i-1)^2/rho*tl_h(i-1);
  ad_eta(i-1)=ad_eta(i-1)+ 1                           *ad_eta(i);
  ad_Sxx(i-1)=ad_Sxx(i-1)+ 1/g/h(i-1)/rho              *ad_eta(i);
  ad_Sxx(i)  =ad_Sxx(i)  - 1/g/h(i-1)/rho              *ad_eta(i);
  ad_h(i-1)  =ad_h(i-1)  - (Sxx(i-1)-Sxx(i))/g/h(i-1)^2/rho*ad_eta(i);
  ad_eta(i)=0;
  %58 tl_Sxx(i) = ...
  %     + tl_Ew(i)*term1 ...
  %     + Ew(i)*tl_term1 ...
  %     + 2*tl_Er(i)*term2 ...
  %     + 2*Er(i)*tl_term2;
  ad_Ew(i) =ad_Ew(i) + term1  *ad_Sxx(i);
  ad_term1=ad_term1+ Ew(i)  *ad_Sxx(i);
  ad_Er(i)=ad_Er(i)+ 2*term2*ad_Sxx(i);
  ad_term2=ad_term2+ 2*Er(i)*ad_Sxx(i);
  ad_Sxx(i)=0;
  %57 tl_term2 = 2*cos(theta(i))*sin(theta(i))*tl_theta(i);
  ad_theta(i)=ad_theta(i)+ 2*cos(theta(i))*sin(theta(i))*ad_term2;
  ad_term2=0;
  %56 tl_term1 = ...
  %     - 2*cos(theta(i))*sin(theta(i))*cg(i)/c(i)*tl_theta(i) ...
  %     + (cos(theta(i))^2+1)*tl_cg(i)/c(i) ...
  %     - (cos(theta(i))^2+1)*cg(i)/c(i)^2*tl_c(i);
  ad_theta(i)=ad_theta(i)- 2*cos(theta(i))*sin(theta(i))*cg(i)/c(i)*ad_term1;
  ad_cg(i)   =ad_cg(i)   + (cos(theta(i))^2+1)/c(i)                *ad_term1;
  ad_c(i)    =ad_c(i)    - (cos(theta(i))^2+1)*cg(i)/c(i)^2        *ad_term1;
  ad_term1=0;

  % update wave height
  %55 tl_H(i) = .5/sqrt(8/rho/g*Ew(i))*8/rho/g*tl_Ew(i);
  ad_Ew(i)=ad_Ew(i)+ .5/sqrt(8/rho/g*Ew(i))*8/rho/g*ad_H(i);
  ad_H(i)=0;
  if(Ew(i)<.001)
    %54 tl_Ew(i)=0;
    ad_Ew(i)=0;
  end

  % update roller energy
  if(~strcmp(betaType,'none'))
    nums=2*Er(i-1)*c(i-1)*cos(theta(i-1))+dx*(Db(i-1)-Dr(i-1));
    denoms=2*c(i)*cos(theta(i));
    if(Er(i)<0)
      %53a1 tl_Er(i)=0;
      ad_Er(i)=0;
    else
      %53b1 tl_Er(i) = tl_nums/denoms - nums/denoms^2*tl_denoms;
      ad_nums  =ad_nums  + 1/denoms     *ad_Er(i);
      ad_denoms=ad_denoms- nums/denoms^2*ad_Er(i);
      ad_Er(i)=0;
    end
    %52 tl_denoms = ...
    %     + 2*tl_c(i)*cos(theta(i)) ...
    %     - 2*c(i)*sin(theta(i))*tl_theta(i);
    ad_c(i)    =ad_c(i)    + 2*cos(theta(i))     *ad_denoms;
    ad_theta(i)=ad_theta(i)- 2*c(i)*sin(theta(i))*ad_denoms;
    ad_denoms=0;
    %51 tl_nums = ...
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
    %50 tl_Dr(i-1) = ...
    %     + 2*g*tl_Er(i-1)*sin(beta(i-1))/c(i-1) ...
    %     + 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)*tl_beta(i-1) ...
    %     - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*tl_c(i-1);
    ad_Er(i-1)  =ad_Er(i-1)  + 2*g*sin(beta(i-1))/c(i-1)          *ad_Dr(i-1);
    ad_beta(i-1)=ad_beta(i-1)+ 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)  *ad_Dr(i-1);
    ad_c(i-1)   =ad_c(i-1)   - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*ad_Dr(i-1);
    ad_Dr(i-1)=0;
    if(strcmp(betaType,'rafati21'))  % rafati et al. (2021) variable-beta
      if(khtot<0.45)
        %49a1 tl_beta(i-1)=0;
        ad_beta(i-1)=0;
      else
        %49b1 tl_beta(i-1) = ...
        %     + 0.03*tl_khtot*(htot(i-1)-H(i-1))/H(i-1) ...
        %     + 0.03*khtot*(tl_htot(i-1)-tl_H(i-1))/H(i-1) ...
        %     - 0.03*khtot*(htot(i-1)-H(i-1))/H(i-1)^2*tl_H(i-1);
        ad_khtot    =ad_khtot    + 0.03*(htot(i-1)-H(i-1))/H(i-1)        *ad_beta(i-1);
        ad_htot(i-1)=ad_htot(i-1)+ 0.03*khtot/H(i-1)                     *ad_beta(i-1);
        ad_H(i-1)   =ad_H(i-1)   - 0.03*khtot/H(i-1)                     *ad_beta(i-1);
        ad_H(i-1)   =ad_H(i-1)   - 0.03*khtot*(htot(i-1)-H(i-1))/H(i-1)^2*ad_beta(i-1);
        ad_beta(i-1)=0;
        if(beta(i-1)>0.1)
          %49b1a1 tl_beta(i-1)=0;
          ad_beta(i-1)=0;
        end
      end
    end
  end  % roller

  % update wave energy
  nums1=cg(i-1)*Ew(i-1)*cos(theta(i-1));
  nums2=Db(i-1)*dx;
  denoms=cg(i-1)*cos(theta(i-1));
  %48 tl_Ew(i) = ...
  %     + (tl_nums1-tl_nums2)/denoms ...
  %     - (nums1-nums2)/denoms^2*tl_denoms;
  ad_nums1 =ad_nums1 + 1/denoms              *ad_Ew(i);
  ad_nums2 =ad_nums2 - 1/denoms              *ad_Ew(i);
  ad_denoms=ad_denoms- (nums1-nums2)/denoms^2*ad_Ew(i);
  ad_Ew(i)=0;
  %47 tl_denoms = ...
  %     + tl_cg(i-1)*cos(theta(i-1)) ...
  %     - cg(i-1)*sin(theta(i-1))*tl_theta(i-1);
  ad_cg(i-1)   =ad_cg(i-1)   + cos(theta(i-1))        *ad_denoms;
  ad_theta(i-1)=ad_theta(i-1)- cg(i-1)*sin(theta(i-1))*ad_denoms;
  ad_denoms=0;
  %46 tl_nums2=tl_Db(i-1)*dx;
  ad_Db(i-1)=ad_Db(i-1)+dx*ad_nums2;
  ad_nums2=0;
  %45 tl_nums1 = ...
  %     + tl_cg(i-1)*Ew(i-1)*cos(theta(i-1)) ...
  %     + cg(i-1)*tl_Ew(i-1)*cos(theta(i-1)) ...
  %     - cg(i-1)*Ew(i-1)*sin(theta(i-1))*tl_theta(i-1);
  ad_cg(i-1)   =ad_cg(i-1)   + Ew(i-1)*cos(theta(i-1))        *ad_nums1;
  ad_Ew(i-1)   =ad_Ew(i-1)   + cg(i-1)*cos(theta(i-1))        *ad_nums1;
  ad_theta(i-1)=ad_theta(i-1)- cg(i-1)*Ew(i-1)*sin(theta(i-1))*ad_nums1;
  ad_nums1=0;

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i-1)/Hm(i-1);
  if(B<=.5)
    Qo=0;
  else
    Qo=(2*B-1)^2;
  end
  c1=alpha/4*rho*g*(omega/2/pi);
  %44 tl_Db(i-1) = ...
  %     + tl_c1*Qb(i-1)*Hm(i-1)^2 ...
  %     + c1*tl_Qb(i-1)*Hm(i-1)^2 ...
  %     + 2*c1*Qb(i-1)*Hm(i-1)*tl_Hm(i-1);
  ad_c1     =ad_c1     + Qb(i-1)*Hm(i-1)^2   *ad_Db(i-1);
  ad_Qb(i-1)=ad_Qb(i-1)+ c1*Hm(i-1)^2        *ad_Db(i-1);
  ad_Hm(i-1)=ad_Hm(i-1)+ 2*c1*Qb(i-1)*Hm(i-1)*ad_Db(i-1);
  ad_Db(i-1)=0;
  %43 tl_c1 = alpha/4*rho*g*tl_omega/2/pi;
  ad_omega=ad_omega+alpha/4*rho*g/2/pi*ad_c1;
  ad_c1=0;
  if(B<=.2)
    %42a1 tl_Qb(i-1)=0;
    ad_Qb(i-1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    nums=Qo-exp(args);
    dens=B^2-exp(args);
    %42b4 tl_Qb(i-1) = ...
    %     + tl_Qo ...
    %     - 2*B*nums/dens*tl_B ...
    %     - B^2*tl_nums/dens ...
    %     + B^2*nums/dens^2*tl_dens;
    ad_Qo  =ad_Qo  + 1              *ad_Qb(i-1);
    ad_B   =ad_B   - 2*B*nums/dens  *ad_Qb(i-1);
    ad_nums=ad_nums- B^2/dens       *ad_Qb(i-1);
    ad_dens=ad_dens+ B^2*nums/dens^2*ad_Qb(i-1);
    ad_Qb(i-1)=0;
    %42b3 tl_dens = 2*B*tl_B - exp(args)*tl_args;
    ad_B   =ad_B   + 2*B      *ad_dens;
    ad_args=ad_args- exp(args)*ad_dens;
    ad_dens=0;
    %42b2 tl_nums = tl_Qo ...
    %           - exp(args)*tl_args;
    ad_Qo  =ad_Qo  + 1        *ad_nums;
    ad_args=ad_args- exp(args)*ad_nums;
    ad_nums=0;
    %42b1 tl_args = ...
    %     + tl_Qo/B^2 ...
    %     - 2*(Qo-1)/B^3*tl_B;
    ad_Qo=ad_Qo+ 1/B^2       *ad_args;
    ad_B =ad_B - 2*(Qo-1)/B^3*ad_args;
    ad_args=0;
  else
    %42c1 tl_Qb(i-1)=0;
    ad_Qb(i-1)=0;
  end
  if(B<=.5)
    Qo=0;
    %41a1 tl_Qo=0;
    ad_Qo=0;
  else
    Qo=(2*B-1)^2;
    %41b1 tl_Qo=2*(2*B-1)*2*tl_B;
    ad_B=ad_B+ 2*(2*B-1)*2*ad_Qo;
    ad_Qo=0;
  end
  %40 tl_B = ...
  %     + tl_H(i-1)/Hm(i-1) ...
  %     - H(i-1)/Hm(i-1)^2*tl_Hm(i-1);
  ad_H(i-1) =ad_H(i-1) + 1/Hm(i-1)       *ad_B;
  ad_Hm(i-1)=ad_Hm(i-1)- H(i-1)/Hm(i-1)^2*ad_B;
  ad_B=0;

  % max wave height
  tharg=gamma(i-1)/0.88.*khtot;
  %39 tl_Hm(i-1) = 0.88*( -1/k(i-1)^2*tanh(tharg)*tl_k(i-1) ...
  %                     + 1/k(i-1)*sech(tharg)^2*tl_tharg );
  ad_k(i-1)=ad_k(i-1)- 0.88*1/k(i-1)^2*tanh(tharg)*ad_Hm(i-1);
  ad_tharg =ad_tharg + 0.88/k(i-1)*sech(tharg)^2  *ad_Hm(i-1);
  ad_Hm(i-1)=0;
  %38 tl_tharg=gamma(i-1)/0.88*tl_khtot ...
  %          + tl_gamma(i-1)/0.88*khtot;
  ad_khtot=ad_khtot+ gamma(i-1)/0.88*ad_tharg;
  ad_gamma(i-1)=ad_gamma(i-1)+ 1/0.88*khtot   *ad_tharg;
  ad_tharg=0;

  % gamma can be either calculated based on deep water wave steepness (s0)
  % following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
  % or based on the empirical fit obtained for duck94 by Ruessink et
  % al. (2003).
  %37 tl_gamma(i-1) = tl_gamma(i-1) + tl_dgamma(i-1);
  ad_gamma(i-1) =ad_gamma(i-1) + ad_gamma(i-1);
  ad_dgamma(i-1)=ad_dgamma(i-1)+ ad_gamma(i-1);
  ad_gamma(i-1)=0;
  if(gammaType==2003)
    %36a1 tl_gamma(i-1) = 0.76*tl_khtot;
    ad_khtot=ad_khtot+0.76*ad_gamma(i-1);
    ad_gamma(i-1)=0;
  end

  % define khtot for convenience
  % tl_khtot = ...
  %     + tl_k(i-1)*htot(i-1) ...
  %     + k(i-1)*tl_htot(i-1);
  ad_k(i-1)   =ad_k(i-1)   + htot(i-1)*ad_khtot;
  ad_htot(i-1)=ad_htot(i-1)+ k(i-1)   *ad_khtot;
  ad_khtot=0;

end

% wave refraction, calculated without setup
% tl_theta=1./sqrt(1-(c.*refconst).^2).*( refconst.*tl_c + tl_refconst.*c );
ad_c       =ad_c       + 1./sqrt(1-(c.*refconst).^2).*refconst.*ad_theta;
ad_refconst=ad_refconst+ 1./sqrt(1-(c.*refconst).^2).*c       .*ad_theta;
ad_theta=0;
% tl_refconst = ...
%     + cos(theta0)/c(1)*tl_theta0 ...
%     - sin(theta0)/c(1)^2*tl_c(1);
ad_theta0=ad_theta0+ cos(theta0)/c(1)  *sum(ad_refconst);
ad_c(1)  =ad_c(1)  - sin(theta0)/c(1)^2*sum(ad_refconst);
ad_refconst=0;

% init variables at offshore boundary
term1=(cos(theta(1))^2+1)*cg(1)/c(1)-.5;
term2=cos(theta(1))^2;
%30 tl_refconst = ...
%     + cos(theta0)/c(1)*tl_theta0 ...
%     - sin(theta0)/c(1)^2*tl_c(1);
ad_theta0=ad_theta0+ cos(theta0)/c(1)  *ad_refconst;
ad_c(1)  =ad_c(1)  - sin(theta0)/c(1)^2*ad_refconst;
ad_refconst=0;
%29 tl_htot(1) = tl_h(1) + tl_eta(1);
ad_h(1)  =ad_h(1)  + ad_htot(1);
ad_eta(1)=ad_eta(1)+ ad_htot(1);
ad_htot(1)=0;
%28 tl_eta(1)=0;
ad_eta(1)=0;
%27 tl_Sxx(1) = ...
%     + tl_Ew(1)*term1 ...
%     + Ew(1)*tl_term1 ...
%     + 2*tl_Er(1)*term2 ...
%     + 2*Er(1)*tl_term2;
ad_Ew(1) =ad_Ew(1) + term1  *ad_Sxx(1);
ad_term1=ad_term1+ Ew(1)  *ad_Sxx(1);
ad_Er(1)=ad_Er(1)+ 2*term2*ad_Sxx(1);
ad_term2=ad_term2+ 2*Er(1)*ad_Sxx(1);
ad_Sxx(1)=0;
%26 tl_term2 = -2*cos(theta(1))*sin(theta(1))*tl_theta(1);
ad_theta(1)=ad_theta(1)- 2*cos(theta(1))*sin(theta(1))*ad_term2;
ad_term2=0;
%25 tl_term1 = ...
%     - 2*cos(theta(1))*sin(theta(1))*tl_theta(1)*cg(1)/c(1) ...
%     + (cos(theta(1))^2+1)*tl_cg(1)/c(1) ...
%     - (cos(theta(1))^2+1)*cg(1)/c(1)^2*tl_c(1);
ad_theta(1)=ad_theta(1)- 2*cos(theta(1))*sin(theta(1))*cg(1)/c(1)*ad_term1;
ad_cg(1)   =ad_cg(1)   + (cos(theta(1))^2+1)/c(1)                *ad_term1;
ad_c(1)    =ad_c(1)    - (cos(theta(1))^2+1)*cg(1)/c(1)^2        *ad_term1;
ad_term1=0;
%24 tl_theta(1)=tl_theta0;
ad_theta0 = ad_theta0 + ad_theta(1);
ad_theta(1)=0;
%23 tl_H(1)=tl_H0;
ad_H0=ad_H0+ad_H(1);
ad_H(1)=0;
%22 tl_Er(1)=0;
ad_Er(1)=0;
%21 tl_Ew(1)=rho*g/8*2*H0*tl_H0;
ad_H0=ad_H0+rho*g/8*2*H0*ad_Ew(1);
ad_Ew(1)=0;

% wave dispersion, calculated without consideration of wave-induced setup
kh=k.*h;
ad_kh=0;  % init AD
%35 tl_cg = tl_n.*c + n.*tl_c;
ad_n=ad_n+ c.*ad_cg;
ad_c=ad_c+ n.*ad_cg;
ad_cg=0;
%34 tl_n = tl_kh./sinh(2*kh) ...
%     - kh./sinh(2*kh).^2.*cosh(2*kh)*2.*tl_kh;
ad_kh=ad_kh+ 1./sinh(2*kh)                      .*ad_n;
ad_kh=ad_kh- kh./sinh(2*kh).^2.*cosh(2*kh)*2.*ad_n;
ad_n=0;
%33 tl_c = ...
%     + tl_omega./k ...
%     - omega./k.^2.*tl_k;
ad_omega =ad_omega + sum(1./k.*ad_c);
ad_k=ad_k- omega./k.^2.*ad_c;
ad_c=0;
%32 tl_k = ...
%     - tl_h.*k.^2.*sech(kh).^2./(tanh(kh)+kh.*sech(kh).^2) ...
%     + 2*omega/g./(tanh(kh)+kh.*sech(kh).^2).*tl_omega;
ad_h=ad_h- k.^2.*sech(kh).^2./(tanh(kh)+kh.*sech(kh).^2).*ad_k;
ad_omega    =ad_omega + sum(2*omega/g./(tanh(kh)+kh.*sech(kh).^2).*ad_k);
ad_k=0;

% init all variables for explicit stepping scheme
%20 tl_Qb   =zeros(nx,1);
%15 tl_theta=zeros(nx,1);
%14 tl_htot =zeros(nx,1);
%13 tl_eta  =zeros(nx,1);
%12 tl_H    =zeros(nx,1);
%12 tl_Dr   =zeros(nx,1);
%11 tl_Db   =zeros(nx,1);
%10 tl_Er   =zeros(nx,1);
%9 tl_Ew   =zeros(nx,1);
ad_Qb   =zeros(nx,1);
ad_theta=zeros(nx,1);
ad_htot =zeros(nx,1);
ad_eta  =zeros(nx,1);
ad_Sxx  =zeros(nx,1);
ad_Hm   =zeros(nx,1);
ad_H    =zeros(nx,1);
ad_Dr   =zeros(nx,1);
ad_Db   =zeros(nx,1);
ad_Er   =zeros(nx,1);
ad_Ew   =zeros(nx,1);

if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);
  s0=H0/L0;
  %8a1 tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0*ones(nx,1);
  ad_s0=ad_s0+ 0.4*sech(33*s0).^2.*33*nx;
else
  %8b1 tl_gamma=zeros(nx,1);  % init
  ad_gamma=zeros(nx,1);  % init
end
if(~exist('dgamma'))
  %7a1 tl_dgamma=zeros(nx,1);  % for compatibility
  ad_dgamma=zeros(nx,1);  % for compatibility
end
if(strcmp(betaType,'none'))
  %6a1 tl_Dr=zeros(nx,1);
  ad_Dr=zeros(nx,1);
end
if(strcmp(betaType,'const'))
  %5a1 tl_beta=zeros(nx,1);
  ad_beta=zeros(nx,1);
end
%4 tl_beta=zeros(nx,1);  % init
ad_beta=zeros(nx,1);  % init

% tl_tauw=tl_tauwin(:,2);  % for compatibility
ad_tauwin(:,2)=ad_tauwin(:,2)+ad_tauw;
ad_tauw=0;

