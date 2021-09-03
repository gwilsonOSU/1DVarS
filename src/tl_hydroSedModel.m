function [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hp] = ...
    tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,...
                     tl_dgamma,...
                     tl_tau_wind,tl_detady,tl_d50,tl_d90,tl_params,bkgd)
%
% TL-code for hydroSedModel.m
%
% note, 'bkgd' struct should be taken from last output of
% hydroSedModel.m, this contains all the NL calculated variables
%

physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

nx=length(h);

%----------------------------------------
% begin TL code
%----------------------------------------

% min depth constraint
h(imask)=hmin;
Q(imask)=0;
Qx(imask)=0;
tl_h(imask)=0;

% % wind stress, following Reniers et al. (2004) eqn (6)
% cd=bkgd.cd;
% w1=windW(:,1).^2+windW(:,2).^2;
% tl_w1=2*windW(:,1).*tl_windW(:,1) ...
%       + 2*windW(:,2).*tl_windW(:,2);
% ind=find(w1>0);
% tl_tau_wind=zeros(nx,2);
% tl_tau_wind(ind,1)=cd*rhoa*( .5./sqrt(w1(ind)).*windW(ind,1).*tl_w1(ind) ...
%                         + sqrt(w1(ind)).*tl_windW(ind,1) );
% tl_tau_wind(ind,2)=cd*rhoa*( .5./sqrt(w1(ind)).*windW(ind,2).*tl_w1(ind) ...
%                         + sqrt(w1(ind)).*tl_windW(ind,2) );

% 1DH wave and longshore current balance
% [Hrms,theta,vbar,kabs,Ew,Er,Dr,hydro_bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma);
[tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
    tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,tl_detady,tl_dgamma,...
                  bkgd.hydro_bkgd);

% convert from (kabs,theta) to vector wavenumber, and calculate c for
% convenience
tl_c = -omega./kabs.^2.*tl_kabs ...
       + tl_omega./kabs;
tl_kvec(:,1)=tl_kabs.*cos(theta)-tl_theta.*kabs.*sin(theta);
tl_kvec(:,2)=tl_kabs.*sin(theta)+tl_theta.*kabs.*cos(theta);

% depth averaged mean flow, Nx2 vector
tl_ubar=-(tl_Ew+2*tl_Er)./(rho*c.*h) ...
        + rho*(Ew+2*Er)./(rho*c.*h).^2.*(tl_c.*h+tl_h.*c);
tl_ubarvec(:,1)=tl_ubar;
tl_ubarvec(:,2)=tl_vbar;

% settling velocity: use Brown & Lawler
tl_ws=tl_ws_brownLawler(tl_d50,d50);

% Reniers et al. (2004) model for velocity at top of boundary layer
if(strcmp(bkgd.sedmodel,'dubarbier') | strcmp(bkgd.sedmodel,'vanderA'))
  tl_udelta=zeros(nx,2);
  tl_delta=zeros(nx,1);
  for i=1:nx
    if(Dr(i)==0)
      tl_udelta(i,:)=[0 0];
      tl_delta(i)=0;
    else
      [tl_udelta(i,:),tl_delta(i)] = ...
          tl_udelta_reniers2004(tl_ubarvec(i,:),tl_kvec(i,:),...
                                tl_omega,tl_h(i),tl_Hrms(i),tl_detady(i),...
                                tl_tau_wind(i,:),tl_Dr(i),tl_params.fv,...
                                tl_d50(i),udel_bkgd(i));
    end
  end
end

% rotate udelta into wave direction, as assumed by sed transport equations
tl_udelta_w(:,1) = ...
    + tl_udelta(:,1).*cos(theta) ...
    - udelta(:,1).*sin(theta).*tl_theta ...
    - tl_udelta(:,2).*sin(theta) ...
    - udelta(:,2).*cos(theta).*tl_theta;
tl_udelta_w(:,2) = ...
    + tl_udelta(:,1).*sin(theta) ...
    + udelta(:,1).*cos(theta).*tl_theta ...
    + tl_udelta(:,2).*cos(theta) ...
    - udelta(:,2).*sin(theta).*tl_theta;

% transport model
if(strcmp(bkgd.sedmodel,'dubarbier'))  % Dubarbier et al. (2015)
  tl_tanbeta=tl_calcTanbeta(tl_h,x)';
  tl_Q=tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta_w,tl_ws,...
                           tl_params.Cw,tl_params.Cc,tl_params.Cf,tl_params.Ka,...
                           bkgd_qtrans);
elseif(strcmp(bkgd.sedmodel,'soulsbyVanRijn'))  % Soulsby & van Rijn
  tl_tanbeta=tl_calcTanbeta(tl_h,x)';
  tl_Q=tl_qtrans_soulsbyVanRijn(tl_d50,tl_d90,tl_h,tl_tanbeta,...
                                tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubarvec,...
                                tl_Dr,tl_params,bkgd_qtrans);
elseif(strcmp(bkgd.sedmodel,'vanderA'))  % van Der A et al. (2013)
  tl_Q=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,tl_omega,...
                         tl_udelta_w,tl_delta,tl_ws,tl_params,bkgd_qtrans);
end

% mitigate transport discontinuity at the shoreline
Q(imask)=0;
tl_Q(imask)=0;

% rotate output from wave-following coords to cartesian
tl_Qx = tl_Q.*cos(theta) ...
        - Q.*sin(theta).*tl_theta;

% bathymetry update: dhdt = -dzdt = dQdx.  This is the Exner equation,
% e.g. see Dubarbier et al. (2015) eqn. (16), and note Q is the volumetric
% transport rate (m2/s) including the bed porosity
if(doMarieu)  % use "stable" Marieu formulation for dh/dt
  dx=abs(x(2)-x(1));
  [tl_hp1,tl_qp1]=tl_dhdt_marieu2007(tl_Qx,tl_h,bkgd_marieu_step1);
  [tl_hp,tl_qp]=tl_dhdt_marieu2007(tl_qp1,tl_hp1,bkgd_marieu_step2);
  tl_dh=tl_hp-tl_h;
  tl_dQdx=nan(nx,1);
else
  tl_dQdx = tl_ddx_upwind(tl_Qx,x,Qx,horig);
  tl_dh=tl_dQdx*dt;  % use ddx_upwind() result
  tl_qp=nan(nx,1);
end
tl_dh=tl_dh.*wgt;   % apply damping near shore
tl_dh(isnan(dh))=0;
tl_hp = tl_h + tl_dh;

% % TEST
% eval(['tl_Q=tl_' outvar ';']);
