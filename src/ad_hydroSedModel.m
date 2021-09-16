function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_dgamma,...
    ad_tau_wind,ad_detady,ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hpout,bkgd)

nx=length(bkgd(1).x);

% init AD
ad_H0=0;
ad_theta0=0;
ad_omega=0;
ad_ka_drag=0;
ad_dgamma=zeros(nx,1);
ad_tau_wind=zeros(nx,1);
ad_detady=zeros(nx,1);
ad_d50=zeros(nx,1);
ad_d90=zeros(nx,1);
ad_h=zeros(nx,1);
ad_hp=zeros(nx,bkgd(1).nsubsteps+1);
ad_params.fv =0;
if(strcmp(bkgd(1).sedmodel,'dubarbier'))
  ad_params.Cw   =0;
  ad_params.Cc   =0;
  ad_params.Cf   =0;
  ad_params.Ka   =0;
elseif(strcmp(bkgd(1).sedmodel,'vanderA'))
  ad_params.n    =0;
  ad_params.m    =0;
  ad_params.xi   =0;
  ad_params.alpha=0;
elseif(strcmp(bkgd(1).sedmodel,'soulsbyVanRijn'))
  ad_params.alphab=0;
  ad_params.facua =0;
end

%---------------------------------------
% begin AD
%---------------------------------------

% drop initial condition from hpout
for n=bkgd(1).nsubsteps:-1:1
  % tl_hpout(:,n)=tl_hp(:,n+1);
  ad_hp(:,n+1) = ad_hp(:,n+1) + ad_hpout(:,n);
  ad_hpout(:,n)=0;
end

% sub-stepping loop
for n=bkgd(1).nsubsteps:-1:1

  [ad_hp1,ad_H01,ad_theta01,ad_omega1,ad_ka_drag1,ad_dgamma1,...
   ad_tau_wind1,ad_detady1,ad_d501,ad_d901,ad_params1] = ...
      ad_hydroSedModel_main(ad_Hrms(:,n),ad_vbar(:,n),ad_theta(:,n),...
                            ad_kabs(:,n),ad_Qx(:,n),ad_hp(:,n+1),bkgd(n));

  ad_hp(:,n) = ad_hp(:,n) + ad_hp1;
  ad_H0 = ad_H0 + ad_H01;
  ad_theta0 = ad_theta0 + ad_theta01;
  ad_omega = ad_omega + ad_omega1;
  ad_ka_drag = ad_ka_drag + ad_ka_drag1;
  ad_dgamma = ad_dgamma + ad_dgamma1;
  ad_tau_wind = ad_tau_wind + ad_tau_wind1;
  ad_detady = ad_detady + ad_detady1;
  ad_d50 = ad_d50 + ad_d501;
  ad_d90 = ad_d90 + ad_d901;
  ad_params.fv = ad_params.fv + ad_params1.fv;
  if(strcmp(bkgd(1).sedmodel,'dubarbier'))
    ad_params.Cw   =ad_params.Cw + ad_params1.Cw;
    ad_params.Cc   =ad_params.Cc + ad_params1.Cc;
    ad_params.Cf   =ad_params.Cf + ad_params1.Cf;
    ad_params.Ka   =ad_params.Ka + ad_params1.Ka;
  elseif(strcmp(bkgd(1).sedmodel,'vanderA'))
    ad_params.n    =ad_params.n     + ad_params1.n    ;
    ad_params.m    =ad_params.m     + ad_params1.m    ;
    ad_params.xi   =ad_params.xi    + ad_params1.xi   ;
    ad_params.alpha=ad_params.alpha + ad_params1.alpha;
  elseif(strcmp(bkgd(1).sedmodel,'soulsbyVanRijn'))
    ad_params.alphab=ad_params.alphab+ad_params1.alphab;
    ad_params.facua =ad_params.facua +ad_params1.facua ;
  end

end  % substepping loop
% tl_hp(:,1) = tl_h;  % init t=0
ad_h = ad_h + ad_hp(:,1);
ad_hp(:,1)=0;

end  % end wrapper function (for sub-stepping loop logic)

% begin main function, for single time step
function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_dgamma,...
    ad_tau_wind,ad_detady,ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel_main(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hp,bkgd)

physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

% min depth constraint
h(imask)=hmin;
Q(imask)=0;
Qx(imask)=0;

nx=length(h);

%----------------------------------------
% begin AD code
%----------------------------------------

% init
ad_h=zeros(nx,1);
ad_dh=zeros(nx,1);
ad_dgamma=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_ka_drag=0;
ad_Ew=zeros(nx,1);
ad_Er=zeros(nx,1);
ad_c=zeros(nx,1);
ad_udelta=zeros(nx,2);
ad_udelta_w=zeros(nx,2);
ad_delta=zeros(nx,1);
ad_ubarvec=zeros(nx,2);
ad_ubar=zeros(nx,1);
ad_detady=zeros(nx,1);
ad_tau_wind=zeros(nx,2);
ad_w1=zeros(nx,1);
ad_Dr=zeros(nx,1);
ad_kvec=zeros(nx,2);
ad_tanbeta=zeros(nx,1);
ad_Q=zeros(nx,1);
ad_dQdx=zeros(nx,1);
ad_d50=zeros(nx,1);
ad_d90=zeros(nx,1);
ad_ws=0;
ad_omega=0;
ad_params.fv =0;
if(strcmp(bkgd.sedmodel,'dubarbier'))
  ad_params.Cw   =0;
  ad_params.Cc   =0;
  ad_params.Cf   =0;
  ad_params.Ka   =0;
elseif(strcmp(bkgd.sedmodel,'vanderA'))
  ad_params.n    =0;
  ad_params.m    =0;
  ad_params.xi   =0;
  ad_params.alpha=0;
elseif(strcmp(bkgd.sedmodel,'soulsbyVanRijn'))
  ad_params.alphab=0;
  ad_params.facua =0;
end
if(doMarieu)
  ad_qp=zeros(nx,1);
end

% special case dt=0: in this case we are done since there is no morphology
% upate to compute.  Skip over all the morphology stuff.
if(dt==0)
  ad_Qx=zeros(nx,1);
  ad_hp=zeros(nx,1);
else

% % bathymetry update: dhdt = -dzdt = dQdx.  This is the Exner equation,
% % e.g. see Dubarbier et al. (2015) eqn. (16), and note Q is the volumetric
% % transport rate (m2/s) including the bed porosity
% tl_hp = tl_h + tl_dh;
ad_h = ad_h + ad_hp;
ad_dh = ad_dh + ad_hp;
ad_hp=0;
% tl_dh(isnan(dh))=0;
ad_dh(isnan(dh))=0;
% tl_dh=tl_dh.*wgt;   % apply damping near shore
ad_dh = ad_dh + wgt.*ad_dh;
if(doMarieu)  % use "stable" Marieu formulation for dh/dt
  dx=abs(x(2)-x(1));
  % tl_dh=tl_hp-tl_h;
  ad_hp = ad_hp + ad_dh;
  ad_h = ad_h - ad_dh;
  ad_dh=0;
  % [tl_hp,tl_qp]=tl_dhdt_marieu2007(tl_qp1,tl_hp1,bkgd_marieu_step2);
  [ad_qp1,ad_hp1]=ad_dhdt_marieu2007(ad_hp,ad_qp,bkgd_marieu_step2);
  ad_hp=0;
  ad_qp=0;
  % [tl_hp1,tl_qp1]=tl_dhdt_marieu2007(tl_Qx,tl_h,bkgd_marieu_step1);
  [ad_Qx,ad_h]=ad_dhdt_marieu2007(ad_hp1,ad_qp1,bkgd_marieu_step1);
  ad_hp1=0;
  ad_qp1=0;
else
  % tl_qp=zeros(nx,1);
  ad_qp=zeros(nx,1);
  % tl_dh=tl_dQdx*dt;  % use ddx_upwind() result
  ad_dQdx = ad_dQdx + dt*ad_dh;
  ad_dh=0;
  % tl_dQdx = tl_ddx_upwind(tl_Qx,x,Qx,horig);
  ad_Qx = ad_ddx_upwind(ad_dQdx,x,Qx,horig);
end

% rotate output from wave-following coords to cartesian
% tl_Qx = tl_Q.*cos(theta) ...
%         - Q.*sin(theta).*tl_theta;
ad_Q    =ad_Q    + cos(theta)   .*ad_Qx;
ad_theta=ad_theta- Q.*sin(theta).*ad_Qx;
ad_Qx=0;

% mitigate transport discontinuity at the shoreline
% tl_Q(imask)=0;
ad_Q(imask)=0;

if(strcmp(bkgd.sedmodel,'dubarbier'))
  %b05 Dubarbier et al. (2015) transport model
  %1 tl_Q=tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta_w,tl_ws,...
  %                          tl_params.Cw,tl_params.Cc,tl_params.Cf,tl_params.Ka,...
  %                          bkgd_qtrans);
  [ad1_tanbeta,ad1_h,ad1_Hrms,ad1_kabs,ad1_omega,ad1_udelta_w,ad1_ws,...
   ad1_Cw,ad1_Cc,ad1_Cf,ad1_Ka] = ...
      ad_qtrans_dubarbier(ad_Q,bkgd_qtrans);
  ad_tanbeta=ad_tanbeta+ad1_tanbeta;
  ad_h      =ad_h      +ad1_h      ;
  ad_Hrms   =ad_Hrms   +ad1_Hrms   ;
  ad_kabs   =ad_kabs   +ad1_kabs   ;
  ad_omega  =ad_omega  +ad1_omega  ;
  ad_udelta_w =ad_udelta_w+ad1_udelta_w;
  ad_ws     =ad_ws     +ad1_ws     ;
  ad_params.Cw=ad_params.Cw+ad1_Cw;
  ad_params.Cc=ad_params.Cc+ad1_Cc;
  ad_params.Cf=ad_params.Cf+ad1_Cf;
  ad_params.Ka=ad_params.Ka+ad1_Ka;
  %0 tl_tanbeta=tl_calcTanbeta(tl_h,x);
  ad_h=ad_h+ad_calcTanbeta(ad_tanbeta,x);
  ad_tanbeta=0*ad_tanbeta;
elseif(strcmp(bkgd.sedmodel,'soulsbyVanRijn'))
  %1 tl_Q=tl_qtrans_soulsbyVanRijn(tl_d50(i),tl_d90,tl_h,tl_tanbeta,...
  %                               tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubarvec,...
  %                               tl_Dr,tl_params,bkgd_qtrans);
  [ad1_d50,ad1_d90,ad1_h,ad1_tanbeta,...
   ad1_Hrms,ad1_kabs,ad1_omega,ad1_theta,ad1_ubarvec,...
   ad1_Dr,ad_params] = ad_qtrans_soulsbyVanRijn(ad_Q,bkgd_qtrans);
  ad_d50     =ad_d50     +ad1_d50     ;
  ad_d90     =ad_d90     +ad1_d90     ;
  ad_h       =ad_h       +ad1_h       ;
  ad_tanbeta =ad_tanbeta+ad1_tanbeta;
  ad_Hrms    =ad_Hrms    +ad1_Hrms    ;
  ad_kabs    =ad_kabs    +ad1_kabs    ;
  ad_omega   =ad_omega   +ad1_omega   ;
  ad_theta   =ad_theta   +ad1_theta   ;
  ad_ubarvec =ad_ubarvec +ad1_ubarvec ;
  ad_Dr      =ad_Dr      +ad1_Dr      ;
  %0 tl_tanbeta=tl_calcTanbeta(tl_h,x)';
  ad_h=ad_h+ad_calcTanbeta(ad_tanbeta,x);
  ad_tanbeta=0*ad_tanbeta;
elseif(strcmp(bkgd.sedmodel,'vanderA'))
  % tl_Q=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,tl_omega,...
  %                        tl_udelta_w,tl_ws,tl_params,bkgd_qtrans);
  [ad1_d50,ad1_d90,ad1_h,ad1_Hrms,ad1_kabs,ad1_omega,...
   ad1_udelta_w,ad1_ws,ad1_params] = ...
      ad_qtrans_vanderA(ad_Q,bkgd_qtrans);
  ad_d50   =ad_d50   +ad1_d50   ;
  ad_d90   =ad_d90   +ad1_d90   ;
  ad_ws    =ad_ws    +ad1_ws    ;
  ad_h=ad_h+ad1_h;
  ad_Hrms  =ad_Hrms  +ad1_Hrms  ;
  ad_kabs  =ad_kabs  +ad1_kabs  ;
  ad_omega =ad_omega +ad1_omega ;
  ad_udelta_w=ad_udelta_w+ad1_udelta_w;
  ad_params.n    =ad_params.n    +ad1_params.n    ;
  ad_params.m    =ad_params.m    +ad1_params.m    ;
  ad_params.xi   =ad_params.xi   +ad1_params.xi   ;
  ad_params.alpha=ad_params.alpha+ad1_params.alpha;
  ad_Q=0;
end

% rotate udelta into wave direction, as assumed by sed transport equations
% tl_udelta_w(:,2) = ...
%     + tl_udelta(:,1).*sin(theta) ...
%     + udelta(:,1).*cos(theta).*tl_theta ...
%     + tl_udelta(:,2).*cos(theta) ...
%     - udelta(:,2).*sin(theta).*tl_theta;
ad_udelta(:,1)=ad_udelta(:,1)+ sin(theta)             .*ad_udelta_w(:,2);
ad_theta      =ad_theta      + udelta(:,1).*cos(theta).*ad_udelta_w(:,2);
ad_udelta(:,2)=ad_udelta(:,2)+ cos(theta)             .*ad_udelta_w(:,2);
ad_theta      =ad_theta      - udelta(:,2).*sin(theta).*ad_udelta_w(:,2);
ad_udelta_w(:,2)=0;
% tl_udelta_w(:,1) = ...
%     + tl_udelta(:,1).*cos(theta) ...
%     - udelta(:,1).*sin(theta).*tl_theta ...
%     - tl_udelta(:,2).*sin(theta) ...
%     - udelta(:,2).*cos(theta).*tl_theta;
ad_udelta(:,1)=ad_udelta(:,1)+ cos(theta)             .*ad_udelta_w(:,1);
ad_theta      =ad_theta      - udelta(:,1).*sin(theta).*ad_udelta_w(:,1);
ad_udelta(:,2)=ad_udelta(:,2)- sin(theta)             .*ad_udelta_w(:,1);
ad_theta      =ad_theta      - udelta(:,2).*cos(theta).*ad_udelta_w(:,1);
ad_udelta_w(:,1)=0;

%b04 Reniers et al. (2004) model for velocity at top of boundary layer
if(strcmp(bkgd.sedmodel,'dubarbier') | strcmp(bkgd.sedmodel,'vanderA'))
  for i=nx:-1:1
    if(Dr(i)==0)
      % tl_udelta(i,:)=[0 0];
      ad_udelta(i,:)=0;
      % tl_delta(i)=0;
      ad_delta(i)=0;
    else
      % [tl_udelta(i,:),tl_delta(i)]=tl_udelta_reniers2004(tl_ubarvec(i,:),tl_kvec(i,:),...
      %                                      tl_h(i),tl_Hrms(i),tl_detady(i),...
      %                                      tl_windW(i,:),tl_Dr(i),tl_fv,tl_d50(i),...
      %                                      udel_bkgd(i));
      [ad1_ubarvec,ad1_kvec,ad1_omega,...
       ad1_h,ad1_Hrms,ad1_detady,...
       ad1_tau_wind,ad1_Dr,ad1_fv,ad1_d50(i)] = ...
          ad_udelta_reniers2004(ad_udelta(i,:),ad_delta(i),udel_bkgd(i));
      ad_ubarvec(i,:)=ad_ubarvec(i,:)+ad1_ubarvec;
      ad_kvec(i,:)   =ad_kvec(i,:)   +ad1_kvec;
      ad_omega       =ad_omega       +ad1_omega;
      ad_h(i)        =ad_h(i)        +ad1_h;
      ad_Hrms(i)     =ad_Hrms(i)     +ad1_Hrms;
      ad_detady(i)   =ad_detady(i)   +ad1_detady;
      ad_tau_wind(i,:)=ad_tau_wind(i,:)+ad1_tau_wind;
      ad_Dr(i)       =ad_Dr(i)       +ad1_Dr;
      ad_params.fv   =ad_params.fv   +ad1_fv;
      ad_d50(i)      =ad_d50(i)      +ad1_d50(i);
    end
  end
else
  ad_params.fv=0;
end

%b03 settling velocity: use Brown & Lawler

%1 tl_ws=tl_ws_brownLawler(tl_d50,d50);
ad_d50=ad_d50+ad_ws_brownLawler(ad_ws,d50);
ad_ws=0;

%b02 depth averaged mean flow, Nx2 vector

%3 tl_ubarvec(:,2)=tl_vbar;
ad_vbar=ad_vbar+ad_ubarvec(:,2);
ad_ubarvec(:,2)=0;
%2 tl_ubarvec(:,1)=tl_ubar;
ad_ubar=ad_ubar+ad_ubarvec(:,1);
ad_ubarvec(:,1)=0;

%1 tl_ubar=-(tl_Ew+2*tl_Er)./(rho*c.*h) ...
%         + rho*(Ew+2*Er)./(rho*c.*h).^2.*(tl_c.*h+tl_h.*c);
ad_Ew=ad_Ew- 1./(rho*c.*h)                  .*ad_ubar;
ad_Er=ad_Er- 2./(rho*c.*h)                  .*ad_ubar;
ad_c =ad_c + rho*(Ew+2*Er)./(rho*c.*h).^2.*h.*ad_ubar;
ad_h =ad_h + rho*(Ew+2*Er)./(rho*c.*h).^2.*c.*ad_ubar;
ad_ubar=0;

%b01 convert from (k,theta) to vector wavenumber, and calculate c for
% convenience

%2b tl_kvec(:,2)=tl_kabs.*sin(theta)+tl_theta.*kabs.*cos(theta);
ad_kabs =ad_kabs +sin(theta)      .*ad_kvec(:,2);
ad_theta=ad_theta+kabs.*cos(theta).*ad_kvec(:,2);
%2a tl_kvec(:,1)=tl_kabs.*cos(theta)-tl_theta.*kabs.*sin(theta);
ad_kabs =ad_kabs +cos(theta)      .*ad_kvec(:,1);
ad_theta=ad_theta-kabs.*sin(theta).*ad_kvec(:,1);
ad_kvec=0;
%1 tl_c=-omega./kabs.^2.*tl_kabs ...
%       + tl_omega./kabs;
ad_kabs=ad_kabs-omega./kabs.^2.*ad_c;
ad_omega=ad_omega+ sum(ad_c./kabs);
ad_c=0;

end  % catch for special case dt==0

%b00 1DH wave and longshore current balance

%1 [tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
%     tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,...
%                   tl_tau_wind,tl_detady,bkgd.hydro_bkgd);
[ad1_h,ad1_H0,ad1_theta0,ad1_omega,ad1_ka_drag,ad1_tauw,ad1_detady,ad1_dgamma] = ...
    ad_hydro_ruessink2001(ad_Hrms,ad_theta,ad_vbar,ad_kabs,ad_Ew,ad_Er,ad_Dr,bkgd.hydro_bkgd);
ad_h      =ad_h      +ad1_h      ;
ad_H0     =ad_H0     +ad1_H0     ;
ad_theta0 =ad_theta0 +ad1_theta0 ;
ad_omega  =ad_omega  +ad1_omega  ;
ad_ka_drag=ad_ka_drag+ad1_ka_drag;
ad_tau_wind=ad_tau_wind+ad1_tauw;
ad_detady=ad_detady+ad1_detady;
ad_dgamma=ad_dgamma+ad1_dgamma;
ad_Hrms =0;
ad_theta=0;
ad_vbar =0;
ad_kabs =0;
ad_Ew=0;
ad_Er=0;
ad_Dr=0;

% % wind stress, following Reniers et al. (2004) eqn (6)
% cd=bkgd.cd;
% w1=windW(:,1).^2+windW(:,2).^2;
% ad_windW=zeros(nx,2);
% ind=find(w1>0);
% % tl_tau_wind(ind,2)=cd*rhoa*( .5./sqrt(w1(ind)).*windW(ind,2).*tl_w1(ind) ...
% %                         + sqrt(w1(ind)).*tl_windW(ind,2) );
% ad_w1(ind)=ad_w1(ind) + cd*rhoa*.5./sqrt(w1(ind)).*windW(ind,2).*ad_tau_wind(ind,2);
% ad_windW(ind,2)=ad_windW(ind,2) + cd*rhoa*sqrt(w1(ind)).*ad_tau_wind(ind,2);
% % tl_tau_wind(ind,1)=cd*rhoa*( .5./sqrt(w1(ind)).*windW(ind,1).*tl_w1(ind) ...
% %                         + sqrt(w1(ind)).*tl_windW(ind,1) );
% ad_w1(ind)=ad_w1(ind) + cd*rhoa*.5./sqrt(w1(ind)).*windW(ind,1).*ad_tau_wind(ind,1);
% ad_windW(ind,1)=ad_windW(ind,1) + cd*rhoa*sqrt(w1(ind)).*ad_tau_wind(ind,1);
% ad_tau_wind=0;

% min depth constraint
% tl_h(imask)=0;
ad_h(imask)=0;

end  % main function
