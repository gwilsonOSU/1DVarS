function [tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr]=tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tauwin,tl_detady,tl_dgamma,bkgd)

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

nx=length(x);
dx=diff(x(1:2));

%---------------------------------------------
% begin TL
%---------------------------------------------

tl_tauw=tl_tauwin(:,2);  % for compatibility

tl_beta=zeros(nx,1);  % init
if(strcmp(betaType,'const'))
  tl_beta=zeros(nx,1);
end
if(strcmp(betaType,'none'))
  tl_Dr=zeros(nx,1);
end
if(~exist('dgamma'))
  tl_dgamma=zeros(nx,1);  % for compatibility
end
if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);
  s0=H0/L0;
  tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0*ones(nx,1);
else
  tl_gamma=zeros(nx,1);  % init
end

% init all variables for explicit stepping scheme
tl_Ew   =zeros(nx,1);
tl_Er   =zeros(nx,1);
tl_Db   =zeros(nx,1);
tl_Dr   =zeros(nx,1);
tl_H    =zeros(nx,1);
tl_Hm   =zeros(nx,1);
tl_Sxx  =zeros(nx,1);
tl_eta  =zeros(nx,1);
tl_htot =zeros(nx,1);
tl_theta=zeros(nx,1);
tl_Qb   =zeros(nx,1);

% wave dispersion, calculated without consideration of wave-induced setup
kh=k.*h;
tl_k = ...
    - tl_h.*k.^2.*sech(kh).^2./(tanh(kh)+kh.*sech(kh).^2) ...
    + 2*omega/g./(tanh(kh)+kh.*sech(kh).^2).*tl_omega;
tl_kh = tl_k.*h + k.*tl_h;
tl_c = ...
    + tl_omega./k ...
    - omega./k.^2.*tl_k;
tl_n = tl_kh./sinh(2*kh) ...
    - kh./sinh(2*kh).^2.*cosh(2*kh)*2.*tl_kh;
tl_cg = tl_n.*c + n.*tl_c;

% init variables at offshore boundary
tl_Ew(1)=rho*g/8*2*H0*tl_H0;
tl_Er(1)=0;
tl_H(1)=tl_H0;
tl_theta(1)=tl_theta0;
term1=(cos(theta(1))^2+1)*cg(1)/c(1)-.5;
tl_term1 = ...
    - 2*cos(theta(1))*sin(theta(1))*tl_theta(1)*cg(1)/c(1) ...
    + (cos(theta(1))^2+1)*tl_cg(1)/c(1) ...
    - (cos(theta(1))^2+1)*cg(1)/c(1)^2*tl_c(1);
term2=cos(theta(1))^2;
tl_term2 = -2*cos(theta(1))*sin(theta(1))*tl_theta(1);
tl_Sxx(1) = ...
    + tl_Ew(1)*term1 ...
    + Ew(1)*tl_term1 ...
    + 2*tl_Er(1)*term2 ...
    + 2*Er(1)*tl_term2;
tl_eta(1)=0;
tl_htot(1) = tl_h(1) + tl_eta(1);

% wave refraction, calculated without setup
tl_refconst = ...
    + cos(theta0)/c(1)*tl_theta0 ...
    - sin(theta0)/c(1)^2*tl_c(1);
tl_theta=1./sqrt(1-(c.*refconst).^2).*( refconst.*tl_c + tl_refconst.*c );

% explicit forward stepping for wave propagation
for i=2:nx

  khtot=k(i-1)*htot(i-1);
  tl_khtot = ...
      + tl_k(i-1)*htot(i-1) ...
      + k(i-1)*tl_htot(i-1);

  % gamma can be either calculated based on deep water wave steepness (s0)
  % following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
  % or based on the empirical fit obtained for duck94 by Ruessink et
  % al. (2003).
  if(gammaType==2003)
    tl_gamma(i-1) = 0.76*tl_khtot;
  end
  tl_gamma(i-1) = tl_gamma(i-1) + tl_dgamma(i-1);

  % max wave height
  tharg=gamma(i-1)/0.88.*khtot;
  tl_tharg=gamma(i-1)/0.88*tl_khtot ...
           + tl_gamma(i-1)/0.88*khtot;
  tl_Hm(i-1) = 0.88*( -1/k(i-1)^2*tanh(tharg)*tl_k(i-1) ...
                      + 1/k(i-1)*sech(tharg)^2*tl_tharg );

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i-1)/Hm(i-1);
  tl_B = ...
      + tl_H(i-1)/Hm(i-1) ...
      - H(i-1)/Hm(i-1)^2*tl_Hm(i-1);
  if(B<=.5)
    Qo=0;
    tl_Qo=0;
  else
    Qo=(2*B-1)^2;
    tl_Qo=2*(2*B-1)*2*tl_B;
  end
  if(B<=.2)
    tl_Qb(i-1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    tl_args = ...
        + tl_Qo/B^2 ...
        - 2*(Qo-1)/B^3*tl_B;
    nums=Qo-exp(args);
    tl_nums = tl_Qo ...
              - exp(args)*tl_args;
    dens=B^2-exp(args);
    tl_dens = 2*B*tl_B - exp(args)*tl_args;
    tl_Qb(i-1) = ...
        + tl_Qo ...
        - 2*B*nums/dens*tl_B ...
        - B^2*tl_nums/dens ...
        + B^2*nums/dens^2*tl_dens;
  else
    tl_Qb(i-1)=0;
  end
  c1=alpha/4*rho*g*(omega/2/pi);
  tl_c1 = alpha/4*rho*g*tl_omega/2/pi;
  tl_Db(i-1) = ...
      + tl_c1*Qb(i-1)*Hm(i-1)^2 ...
      + c1*tl_Qb(i-1)*Hm(i-1)^2 ...
      + 2*c1*Qb(i-1)*Hm(i-1)*tl_Hm(i-1);

  % update wave energy
  nums1=cg(i-1)*Ew(i-1)*cos(theta(i-1));
  tl_nums1 = ...
      + tl_cg(i-1)*Ew(i-1)*cos(theta(i-1)) ...
      + cg(i-1)*tl_Ew(i-1)*cos(theta(i-1)) ...
      - cg(i-1)*Ew(i-1)*sin(theta(i-1))*tl_theta(i-1);
  nums2=Db(i-1)*dx;
  tl_nums2=tl_Db(i-1)*dx;
  denoms=cg(i)*cos(theta(i));
  tl_denoms = ...
      + tl_cg(i)*cos(theta(i)) ...
      - cg(i)*sin(theta(i))*tl_theta(i);
  tl_Ew(i) = ...
      + (tl_nums1-tl_nums2)/denoms ...
      - (nums1-nums2)/denoms^2*tl_denoms;

  % update roller energy
  if(~strcmp(betaType,'none'))
    if(strcmp(betaType,'rafati21'))  % rafati et al. (2021) variable-beta
      if(khtot<0.45)
        tl_beta(i-1)=0;
      else
        tl_beta(i-1) = ...
            + 0.03*tl_khtot*(htot(i-1)-H(i-1))/H(i-1) ...
            + 0.03*khtot*(tl_htot(i-1)-tl_H(i-1))/H(i-1) ...
            - 0.03*khtot*(htot(i-1)-H(i-1))/H(i-1)^2*tl_H(i-1);
        if(beta(i-1)>0.1)
          tl_beta(i-1)=0;
        end
      end
    end
    tl_Dr(i-1) = ...
        + 2*g*tl_Er(i-1)*sin(beta(i-1))/c(i-1) ...
        + 2*g*Er(i-1)*cos(beta(i-1))/c(i-1)*tl_beta(i-1) ...
        - 2*g*Er(i-1)*sin(beta(i-1))/c(i-1)^2*tl_c(i-1);
    nums=2*Er(i-1)*c(i-1)*cos(theta(i-1))+dx*(Db(i-1)-Dr(i-1));
    tl_nums = ...
        + 2*tl_Er(i-1)*c(i-1)*cos(theta(i-1)) ...
        + 2*Er(i-1)*tl_c(i-1)*cos(theta(i-1)) ...
        - 2*Er(i-1)*c(i-1)*sin(theta(i-1))*tl_theta(i-1) ...
        + dx*(tl_Db(i-1)-tl_Dr(i-1));
    denoms=2*c(i)*cos(theta(i));
    tl_denoms = ...
        + 2*tl_c(i)*cos(theta(i)) ...
        - 2*c(i)*sin(theta(i))*tl_theta(i);
    if(Er(i)<0)
      tl_Er(i)=0;
    else
      tl_Er(i) = tl_nums/denoms - nums/denoms^2*tl_denoms;
    end
  end  % roller

  % update wave height
  if(Ew(i)<.001)
    tl_Ew(i)=0;
  end
  tl_H(i) = .5/sqrt(8/rho/g*Ew(i))*8/rho/g*tl_Ew(i);

  % update wave induced setup
  term1=(cos(theta(i))^2+1)*cg(i)/c(i)-.5;
  tl_term1 = ...
      - 2*cos(theta(i))*sin(theta(i))*cg(i)/c(i)*tl_theta(i) ...
      + (cos(theta(i))^2+1)*tl_cg(i)/c(i) ...
      - (cos(theta(i))^2+1)*cg(i)/c(i)^2*tl_c(i);
  term2=cos(theta(i))^2;
  tl_term2 = 2*cos(theta(i))*sin(theta(i))*tl_theta(i);
  tl_Sxx(i) = ...
      + tl_Ew(i)*term1 ...
      + Ew(i)*tl_term1 ...
      + 2*tl_Er(i)*term2 ...
      + 2*Er(i)*tl_term2;
  tl_eta(i) = ...
      tl_eta(i-1) ...
      + (tl_Sxx(i-1)-tl_Sxx(i))/g/htot(i-1)/rho ...
      - (Sxx(i-1)-Sxx(i))/g/htot(i-1)^2/rho*tl_htot(i-1);
  tl_htot(i) = tl_h(i) + tl_eta(i);

end

% radiation stress gradient
if(~strcmp(betaType,'none'))
  tl_dSxydx = ...
      - cos(theta)./c.*Dr/rho.*tl_theta ...
      + sin(theta)./c.^2.*Dr/rho.*tl_c ...
      - sin(theta)./c/rho.*tl_Dr;
else
  tl_dSxydx = ...
      - cos(theta)./c.*Db/rho.*tl_theta ...
      + sin(theta)./c.^2.*Db/rho.*tl_c ...
      - sin(theta)./c/rho.*tl_Db;
end

% total force = radiation stress gradient + wind stress + pressure gradient,
% m2/s2 units
tl_Fy = ...
    + tl_dSxydx ...
    + tl_tauw/rho ...
    + g*tl_htot.*detady ...
    + g*htot.*tl_detady;

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000).  To get TL model, differentiate the eqn for v (i.e., the
% fsolve() line in waveModel.m) on both sides, then solve for
% tl_v
B=a^2+(v./urms).^2;
dens = -urms.*Cd - v.^2.*Cd./urms./B;
tl_Cd=0.015*(1/3)*(ka_drag./htot).^(-2/3).*(-ka_drag./htot.^2.*tl_htot+tl_ka_drag./htot);
tl_urms=1.416*omega*( tl_H./(4*sinh(k.*htot)) ...
                      -H./(4*sinh(k.*htot).^2).*cosh(k.*htot).*( tl_k.*htot+k.*tl_htot ) ) ...
        + 1.416*H./(4*sinh(k.*htot))*tl_omega;
tl_N = tl_Fy./sqrt(B) + tl_urms.*(v.*Cd-v.^3.*Cd./(urms.^2.*B)) + tl_Cd.*v.*urms;
if(nu==0)
  tl_v=tl_N./dens;
else
  tl_v = inv(diag(dens)+A)*tl_N;  % tl_N = (dens + A) * tl_v
end
