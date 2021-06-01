function tl_qtot = tl_qtrans_soulsbyVanRijn(tl_d50,tl_d90,tl_h,tl_tanbeta,...
                                         tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubar,...
                                         tl_Dr,tl_param,bkgd) %,outvar)
%
% TL-code for qtrans_soulsbyVanRijn.m
%

% break out NL background vars
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
psi=bkgd.psi;

physicalConstants;

% drag coef for currents, from soulsby.  He says "z0 should be set to 6mm"
% z0=.006;
% Cd = (.4./(log(h/z0)-1)).^2;
tl_Cd = -2*.4^2./(log(h/z0)-1).^3.*tl_h./h;

% mean flow mag
% uE2 = ubar(:,1).^2 + ubar(:,2).^2;
tl_uE2 = 2*ubar(:,1).*tl_ubar(:,1) ...
         + 2*ubar(:,2).*tl_ubar(:,2);

% critical velocity, Soulsby eqn 133d
% ucr=0.19*d50^.1*log10(4*h./d90);
tl_ucr = .1*0.19*d50^(.1-1)*log10(4*h./d90).*tl_d50 ...
         + 0.19*d50^.1.*d90./(4*h*log(10)).*( ...
             + 4*tl_h./d90 ...
             - 4*h./d90.^2.*tl_d90 );

% rms wave orbital velocity.  Omits xbeach's finite amplitude "correction",
% and turbulence enhancement
% urms = Hrms/2*omega./sinh(kabs.*h);
tl_urms = tl_Hrms/2*omega./sinh(kabs.*h) ...
          - Hrms/2*omega./sinh(kabs.*h).^2.*cosh(kabs.*h).*( ...
              tl_kabs.*h + kabs.*tl_h ) ...
          + Hrms/2./sinh(kabs.*h)*tl_omega;

% wave asymmetry & skewness factors
% Hmo = 1.4*Hrms;
tl_Hmo = 1.4*tl_Hrms;
% aw=Hmo/2;
tl_aw=tl_Hmo/2;
% Ur=3/4*aw.*kabs./(kabs.*h).^3;
tl_Ur = 3/4*tl_aw.*kabs./(kabs.*h).^3 ...
        + 3/4*aw.*tl_kabs./(kabs.*h).^3 ...
        - 3*3/4*aw.*kabs./(kabs.*h).^4.*( ...
            tl_kabs.*h + kabs.*tl_h );
% psi = -pi/2+pi/2*tanh(.64./Ur.^(.60));
tl_psi = -pi/2*sech(.64./Ur.^(.60)).^2.*.6*.64.*Ur.^(-.6-1).*tl_Ur;
% earg = (-.61-log(Ur))/-.35;
tl_earg = 1./Ur/.35.*tl_Ur;
% dens = 1+exp(earg);
tl_dens = exp(earg).*tl_earg;
% Sk = .79./dens.*cos(psi);
tl_Sk = -.79./dens.^2.*cos(psi).*tl_dens ...
        - .79./dens.*sin(psi).*tl_psi;
% As = .79./dens.*sin(psi);
tl_As = -.79./dens.^2.*sin(psi).*tl_dens ...
        +.79./dens.*cos(psi).*tl_psi;

% advection velocity, plus extra term to represent wave asymmetry
% VW = param.facua*urms.*(Sk-As);  % eqn 2.67
tl_VW = tl_param.facua*urms.*(Sk-As) ...
        + param.facua*tl_urms.*(Sk-As) ...
        + param.facua*urms.*(tl_Sk-tl_As);
% uAV = VW.*cos(theta)+ubar(:,1);  % eqn 2.66
tl_uAV = tl_VW.*cos(theta) ...
         - VW.*sin(theta).*tl_theta ...
         + tl_ubar(:,1);

% sediment flux in m2/s, eqns from Soulsby's text
% Dstar=(g*(s-1)/nu^2)^(1/3)*d50;  % SCALAR
tl_Dstar=(g*(s-1)/nu^2)^(1/3)*tl_d50;
% Asb = .005*h.*(d50./h).^1.2/((s-1)*g*d50)^1.2;  % eqn 136b
tl_Asb = .005*tl_h.*(d50./h).^1.2/((s-1)*g*d50)^1.2 ...
         + 1.2*.005*h.*(d50./h).^(1.2-1)/((s-1)*g*d50)^1.2.*( ...
             tl_d50./h - d50./h.^2.*tl_h ) ...
         - 1.2*.005*h.*(d50./h).^1.2.*((s-1)*g*d50)^(-2.2).*(s-1)*g*tl_d50;
% Ass = .012*d50*Dstar^(-.6)/((s-1)*g*d50)^1.2;  % eqn 136c, SCALAR
tl_Ass = .012*tl_d50*Dstar^(-.6)/((s-1)*g*d50)^1.2 ...
         - .6*.012*d50*Dstar^(-.6-1)/((s-1)*g*d50)^1.2.*tl_Dstar ...
         - 1.2*.012*d50*Dstar^(-.6)*((s-1)*g*d50)^(-2.2)*(s-1)*g*tl_d50;
% Ufact = sqrt(uE2 + .018*urms.^2./Cd) - ucr;
tl_Ufact = .5./sqrt(uE2 + .018*urms.^2./Cd).*( ...
    tl_uE2 + 2*.018*urms./Cd.*tl_urms - .018*urms.^2./Cd.^2.*tl_Cd ) ...
    - tl_ucr;
% slopeFact = 1-param.alphab*tanbeta;
tl_slopeFact = -tl_param.alphab*tanbeta ...
    - param.alphab*tl_tanbeta;
% q=zeros(length(x),1);
tl_q=zeros(length(x),1);
% C=zeros(length(x),1);
tl_C=zeros(length(x),1);
ind=find(Ufact>0);
% q(ind) = (Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind);
tl_q(ind) = (tl_Ass+tl_Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind) ...
    + (Ass+Asb(ind)).*tl_uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind) ...
    + 2.4*(Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4-1).*slopeFact(ind).*tl_Ufact(ind) ...
    + (Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*tl_slopeFact(ind);
% C(ind) = (Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind);
tl_C(ind) = (tl_Ass+tl_Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind) ...
    + 2.4*(Ass+Asb(ind)).*Ufact(ind).^(1.4).*slopeFact(ind)./h(ind).*tl_Ufact(ind) ...
    + (Ass+Asb(ind)).*Ufact(ind).^(2.4).*tl_slopeFact(ind)./h(ind) ...
    - (Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind).^2.*tl_h(ind);

% sediment diffusion coefficient, this is copied from xbeach code since it's
% not specified in the manual nor the Roelvink paper
% facDc=1;   % xbeach "turn on sediment diffusion" flag
% nuhfac=1;  % xbeach "turn on roller-induced viscosity" flag
% nuh=.1;    % bkgd viscosity, default 0.1 m2/s in xbeach
% Dh = facDc*(nuh+nuhfac*h.*(Dr/rho).^(1/3));  % copied from xbeach code
ind=find(Dr>0);
tl_Dh=zeros(length(x),1);
tl_Dh(ind) = facDc*nuhfac*tl_h(ind).*(Dr(ind)/rho).^(1/3) ...
        + (1/3)*facDc*nuhfac*h(ind).*(Dr(ind)/rho).^(1/3-1).*tl_Dr(ind)/rho;

% contribution from sediment diffusion
if(length(h)>2)
  % dCdx=ddx_upwind(x,C);
  tl_dCdx=tl_ddx_upwind(tl_C,x,C);
  % qtot = q + Dh.*dCdx;
  tl_qtot = tl_q ...
            + tl_Dh.*ddx_upwind(x,C) ...
            + Dh.*tl_dCdx;
else
  tl_qtot=tl_q;
end

% % TEST: redefine output var
% eval(['tl_qtot = tl_' outvar ';'])
