function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
          ad_detady,ad_dgamma,ad_dAw,ad_dSw,...
          ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hpout,bkgd)%,invar)

nx=length(bkgd(1).x);

% init AD
ad_h=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_omega=0;
ad_ka_drag=0;
ad_beta0=0;
ad_tau_wind=zeros(nx,1);
ad_detady=zeros(nx,1);
ad_dgamma=zeros(nx,1);
ad_dAw=zeros(nx,1);
ad_dSw=zeros(nx,1);
ad_d50=zeros(nx,1);
ad_d90=zeros(nx,1);
ad_params.fv =0;
ad_params.ks =0;
ad_params.lambda=0;
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
  if(isfield(bkgd(1).params,'Cc'))
    ad_params.Cc   =0;
    ad_params.Cf   =0;
  end
elseif(strcmp(bkgd(1).sedmodel,'soulsbyVanRijn'))
  ad_params.alphab=0;
  ad_params.facua =0;
end
ad_hp=zeros(nx,bkgd(1).nsubsteps+1);

%---------------------------------------
% begin AD
%---------------------------------------

% drop initial condition from hpout
for n=bkgd(1).nsubsteps:-1:1
  % tl_hpout(:,n)=tl_hp(:,n+1);
  ad_hp(:,n+1)=ad_hp(:,n+1)+ad_hpout(:,n);
  ad_hpout(:,n)=0;
end

% sub-stepping loop
for n=bkgd(1).nsubsteps:-1:1

  % [tl_Hrms(:,n),tl_vbar(:,n),tl_theta(:,n),tl_kabs(:,n),tl_Qx(:,n),tl_hp(:,n+1)] = ...
  %     tl_hydroSedModel_main(tl_hp(:,n),tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,...
  %                           tl_tau_wind,tl_detady,...
  %                           tl_dgamma,tl_d50,tl_d90,tl_params,...
  %                           bkgd(n),outvar);
  [ad_hp(:,n),ad1_H0,ad1_theta0,ad1_omega,ad1_ka_drag,ad1_beta0,ad1_tau_wind,...
   ad1_detady,ad1_dgamma,ad1_dAw,ad1_dSw,ad1_d50,ad1_d90,ad1_params] = ...
      ad_hydroSedModel_main(ad_Hrms(:,n),ad_vbar(:,n),ad_theta(:,n),...
                            ad_kabs(:,n),ad_Qx(:,n),ad_hp(:,n+1),...
                            bkgd(n));%,invar);

  ad_H0      =ad_H0      +ad1_H0      ;
  ad_theta0  =ad_theta0  +ad1_theta0  ;
  ad_omega   =ad_omega   +ad1_omega   ;
  ad_ka_drag =ad_ka_drag +ad1_ka_drag ;
  ad_beta0   =ad_beta0   +ad1_beta0   ;
  ad_tau_wind=ad_tau_wind+ad1_tau_wind;
  ad_detady  =ad_detady  +ad1_detady  ;
  ad_dgamma  =ad_dgamma  +ad1_dgamma  ;
  ad_dAw     =ad_dAw     +ad1_dAw     ;
  ad_dSw     =ad_dSw     +ad1_dSw     ;
  ad_d50     =ad_d50     +ad1_d50     ;
  ad_d90     =ad_d90     +ad1_d90     ;

  ad_Hrms(:,n) =0;
  ad_vbar(:,n) =0;
  ad_theta(:,n)=0;
  ad_kabs(:,n) =0;
  ad_Qx(:,n)   =0;
  ad_hp(:,n+1) =0;

  ad_params.fv = ad_params.fv + ad1_params.fv;
  ad_params.ks = ad_params.ks + ad1_params.ks;
  ad_params.lambda = ad_params.lambda + ad1_params.lambda;
  if(strcmp(bkgd(1).sedmodel,'dubarbier'))
    ad_params.Cw   =ad_params.Cw + ad1_params.Cw;
    ad_params.Cc   =ad_params.Cc + ad1_params.Cc;
    ad_params.Cf   =ad_params.Cf + ad1_params.Cf;
    ad_params.Ka   =ad_params.Ka + ad1_params.Ka;
  elseif(strcmp(bkgd(1).sedmodel,'vanderA'))
    ad_params.n    =ad_params.n     + ad1_params.n    ;
    ad_params.m    =ad_params.m     + ad1_params.m    ;
    ad_params.xi   =ad_params.xi    + ad1_params.xi   ;
    ad_params.alpha=ad_params.alpha + ad1_params.alpha;
    if(isfield(bkgd(1).params,'Cc'))
      ad_params.Cc   =ad_params.Cc    + ad1_params.Cc   ;
      ad_params.Cf   =ad_params.Cf    + ad1_params.Cf   ;
    end
  elseif(strcmp(bkgd(1).sedmodel,'soulsbyVanRijn'))
    ad_params.alphab=ad_params.alphab+ad1_params.alphab;
    ad_params.facua =ad_params.facua +ad1_params.facua ;
  end

end
% tl_hp(:,1) = tl_h;  % init t=0
ad_h=ad_h+ad_hp(:,1);

end  % end wrapper function (for sub-stepping loop logic)

% begin main function, for single time step
function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
          ad_detady,ad_dgamma,ad_dAw,ad_dSw,ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel_main(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hp,bkgd)%,invar)

physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

nx=length(h);

%----------------------------------------
% init AD
%----------------------------------------

ad_h=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_omega=0;
ad_ka_drag=0;
ad_beta0  =0;
ad_tau_wind=zeros(nx,2);
ad_detady=zeros(nx,1);
ad_dgamma=zeros(nx,1);
ad_dAw=zeros(nx,1);
ad_dSw=zeros(nx,1);
ad_Aw=zeros(nx,1);
ad_Sw=zeros(nx,1);
ad_Aw0=zeros(nx,1);
ad_Sw0=zeros(nx,1);
ad_Uw=zeros(nx,1);
ad_d50=zeros(nx,1);
ad_d90=zeros(nx,1);
ad_params.fv =0;
ad_params.ks =0;
ad_params.lambda =0;
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
  if(isfield(bkgd.params,'Cc'))
    ad_params.Cc   =0;
    ad_params.Cf   =0;
  end
elseif(strcmp(bkgd.sedmodel,'soulsbyVanRijn'))
  ad_params.alphab=0;
  ad_params.facua =0;
end

ad_k=zeros(nx,2);
ad_dh=zeros(nx,1);
ad_Ew=zeros(nx,1);
ad_Er=zeros(nx,1);
ad_Hmo=zeros(nx,1);
ad_c=zeros(nx,1);
ad_udelta=zeros(nx,2);
ad_udelta_w=zeros(nx,2);
ad_delta_bl=zeros(nx,1);
ad_ubarvec=zeros(nx,2);
ad_ubar=zeros(nx,2);
ad_ubar0=zeros(nx,2);
ad_ubarx=zeros(nx,1);
ad_w1=zeros(nx,1);
ad_Dr=zeros(nx,1);
ad_tanbeta=zeros(nx,1);
ad_Q=zeros(nx,1);
ad_Q0=zeros(nx,1);
ad_Q1=zeros(nx,1);
ad_dQdx=zeros(nx,1);
ad_ws=0;
if(doMarieu)
  ad_qp=zeros(nx,1);
end

%----------------------------------------
% begin AD code
%----------------------------------------

% % TEST: look at a specific variable
% ad_vbar =zeros(nx,1);
% ad_theta=zeros(nx,1);
% ad_kabs =zeros(nx,1);
% ad_Hrms =zeros(nx,1);
% ad_Qx   =zeros(nx,1);
% if((length(invar)>=6 & strcmp(invar(1:6),'udelta')) | strcmp(invar,'ubar'))
%   % eval(['tl_hp = tl_' outvar '(:,2);']);
%   eval(['ad_' invar '(:,2)=ad_hp;']);
%   % eval(['tl_' outvar '(:,1)=0;']);
%   eval(['ad_' invar '(:,1)=0;']);
% else
%   % eval(['tl_hp = tl_' outvar ';']);
%   eval(['ad_' invar '=ad_hp;']);
% end
% if(~strcmp(invar,'hp'))
%   ad_hp=zeros(nx,1);
% end

% bathymetry update: dhdt = -dzdt = dQdx.  This is the Exner equation,
% e.g. see Dubarbier et al. (2015) eqn. (16), and note Q is the volumetric
% transport rate (m2/s) including the bed porosity
%21 tl_hp = tl_h + tl_dh;
ad_h =ad_h + ad_hp;
ad_dh=ad_dh+ ad_hp;
ad_hp=0;
%20 tl_dh(isnan(dh))=0;
ad_dh(isnan(dh))=0;
%19 tl_dh = tl_dh.*wgt;   % apply damping near shore
ad_dh=wgt.*ad_dh;
if(doMarieu)  % use "stable" Marieu formulation for dh/dt
  dx=abs(x(2)-x(1));
  %18a4 tl_dQdx=zeros(nx,1);
  ad_dQdx=zeros(nx,1);
  %18a3 tl_dh = tl_hp - tl_h;
  ad_hp=ad_hp+ ad_dh;
  ad_h =ad_h - ad_dh;
  ad_dh=0;
  %18a2 [tl_hp ,tl_qp ]=tl_dhdt_marieu2007(tl_qp1,tl_hp1,bkgd_marieu_step2);
  [ad_qp1,ad_hp1]=ad_dhdt_marieu2007(ad_hp,ad_qp,bkgd_marieu_step2);
  ad_hp=0;
  ad_qp=0;
  %18a1 [tl_hp1,tl_qp1]=tl_dhdt_marieu2007(tl_Qx ,tl_h  ,bkgd_marieu_step1);
  [ad1_Qx,ad1_h]=ad_dhdt_marieu2007(ad_hp1,ad_qp1,bkgd_marieu_step1);
  ad_Qx = ad_Qx + ad1_Qx;
  ad_h = ad_h + ad1_h;
  ad_hp1=0;
  ad_qp1=0;
else
  %18b3 tl_qp = zeros(nx,1);
  ad_qp = zeros(nx,1);
  %18b2 tl_dh = tl_dQdx*dt;
  ad_dQdx=ad_dQdx+ dt*ad_dh;
  ad_dh=0;
  %18b1 tl_dQdx = tl_ddx_upwind(tl_Qx,x,Qx,horig);
  ad_Qx = ad_Qx + ad_ddx_upwind(ad_dQdx,x,Qx,horig);
  ad_dQdx=0;
end

% rotate output from wave-following coords to cartesian
%17 tl_Qx = ...
%     + tl_Q.*cos(theta) ...
%     - Q.*sin(theta).*tl_theta;
ad_Q    =ad_Q    + cos(theta)   .*ad_Qx;
ad_theta=ad_theta- Q.*sin(theta).*ad_Qx;
ad_Qx=0;

% filtering is applied in NL model, neglected in TL
% tl_Q = tl_Q1;
ad_Q1 = ad_Q1 + ad_Q;
ad_Q=0;

% mitigate transport discontinuity at the shoreline
% tl_Q1(imax:end)=0;
ad_Q1(imax:end)=0;
% tl_Q1 = tl_Q0;
ad_Q0 = ad_Q0 + ad_Q1;
ad_Q1=0;

% run the requested model for sediment flux (m2/s)
if(strcmp(sedmodel,'dubarbier'))  % Dubarbier et al. (2015)
  %15a1 tl_Q0 =tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,...
  %                           tl_params.Cw,tl_params.Cc,tl_params.Cf,tl_params.Ka,...
  %                           bkgd_qtrans);
  [ad1_tanbeta,ad1_h,ad1_Hrms,ad1_kabs,ad1_omega,ad1_udelta,ad1_ws,ad1_Aw,ad1_Sw,ad1_Uw,...
   ad1_params_Cw,ad1_params_Cc,ad1_params_Cf,ad1_params_Ka] = ...
      ad_qtrans_dubarbier(ad_Q0,bkgd_qtrans);
  ad_tanbeta  =ad_tanbeta  +ad1_tanbeta  ;
  ad_h        =ad_h        +ad1_h        ;
  ad_Hrms     =ad_Hrms     +ad1_Hrms     ;
  ad_kabs     =ad_kabs     +ad1_kabs     ;
  ad_omega    =ad_omega    +ad1_omega    ;
  ad_udelta   =ad_udelta   +ad1_udelta   ;
  ad_ws       =ad_ws       +ad1_ws       ;
  ad_Aw       =ad_Aw       +ad1_Aw       ;
  ad_Sw       =ad_Sw       +ad1_Sw       ;
  ad_Uw       =ad_Uw       +ad1_Uw       ;
  ad_params.Cw=ad_params.Cw+ad1_params_Cw;
  ad_params.Cc=ad_params.Cc+ad1_params_Cc;
  ad_params.Cf=ad_params.Cf+ad1_params_Cf;
  ad_params.Ka=ad_params.Ka+ad1_params_Ka;
  ad_Q0=0;
elseif(strcmp(sedmodel,'soulsbyVanRijn'))  % Soulsby & van Rijn
  %15b1 tl_Q0 = tl_qtrans_soulsbyVanRijn(tl_d50,tl_d90,tl_h,tl_tanbeta,...
  %                                 tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubar,...
  %                                 tl_Dr,tl_param,bkgd_qtrans);
  [ad1_d50,ad1_d90,ad1_h,ad1_tanbeta,...
   ad1_Hrms,ad1_kabs,ad1_omega,ad1_theta,ad1_ubar,...
   ad1_Dr,ad1_params] = ...
      ad_qtrans_soulsbyVanRijn(ad_Q0,bkgd_qtrans);
  ad_d50         =ad_d50         +ad1_d50         ;
  ad_d90         =ad_d90         +ad1_d90         ;
  ad_h           =ad_h           +ad1_h           ;
  ad_tanbeta     =ad_tanbeta     +ad1_tanbeta     ;
  ad_Hrms        =ad_Hrms        +ad1_Hrms        ;
  ad_kabs        =ad_kabs        +ad1_kabs        ;
  ad_omega       =ad_omega       +ad1_omega       ;
  ad_theta       =ad_theta       +ad1_theta       ;
  ad_ubar        =ad_ubar        +ad1_ubar        ;
  ad_Dr          =ad_Dr          +ad1_Dr          ;
  ad_params.facua =ad_params.facua +ad1_params.facua ;
  ad_params.alphab=ad_params.alphab+ad1_params.alphab;
  ad_Q0=0;
elseif(strcmp(sedmodel,'vanderA'))  % van Der A et al. (2013)
  %15c1 tl_Q0 = tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,...
  %                          tl_udelta,tl_ws,tl_params,bkgd_qtrans);
  [ad1_d50,ad1_d90,ad1_h,ad1_tanbeta,ad1_Hrms,ad1_kabs,ad1_omega,ad1_udelta,ad1_delta_bl,ad1_ws,ad1_Aw,ad1_Sw,ad1_Uw,ad1_params] = ...
      ad_qtrans_vanderA(ad_Q0,bkgd_qtrans);
  ad_d50        =ad_d50        +ad1_d50        ;
  ad_d90        =ad_d90        +ad1_d90        ;
  ad_h          =ad_h          +ad1_h          ;
  ad_tanbeta    =ad_tanbeta    +ad1_tanbeta    ;
  ad_Hrms       =ad_Hrms       +ad1_Hrms       ;
  ad_kabs       =ad_kabs       +ad1_kabs       ;
  ad_omega      =ad_omega      +ad1_omega      ;
  ad_udelta     =ad_udelta     +ad1_udelta     ;
  ad_delta_bl   =ad_delta_bl   +ad1_delta_bl   ;
  ad_ws         =ad_ws         +ad1_ws         ;
  ad_params.n    =ad_params.n    +ad1_params.n    ;
  ad_params.m    =ad_params.m    +ad1_params.m    ;
  ad_params.xi   =ad_params.xi   +ad1_params.xi   ;
  ad_params.alpha=ad_params.alpha+ad1_params.alpha;
  if(isfield(bkgd(1).params,'Cc'))
    ad_params.Cc   =ad_params.Cc   +ad1_params.Cc   ;
    ad_params.Cf   =ad_params.Cf   +ad1_params.Cf   ;
  end
  ad_Aw    =ad_Aw    +ad1_Aw    ;
  ad_Sw    =ad_Sw    +ad1_Sw    ;
  ad_Uw    =ad_Uw    +ad1_Uw    ;
  ad_Q0=0;
end
%14 tl_tanbeta = tl_calcTanbeta(tl_h,x)';
ad_h = ad_h + ad_calcTanbeta(ad_tanbeta,x);
ad_tanbeta=0;

% rotate udelta into wave direction, as assumed by sed transport equations
%13 tl_udelta_w(:,2) = ...
%     - tl_udelta(:,1).*sin(theta) ...
%     - udelta(:,1).*cos(theta).*tl_theta ...
%     + tl_udelta(:,2).*cos(theta) ...
%     - udelta(:,2).*sin(theta).*tl_theta;
ad_udelta(:,1)=ad_udelta(:,1)- sin(theta)             .*ad_udelta_w(:,2);
ad_theta      =ad_theta      - udelta(:,1).*cos(theta).*ad_udelta_w(:,2);
ad_udelta(:,2)=ad_udelta(:,2)+ cos(theta)             .*ad_udelta_w(:,2);
ad_theta      =ad_theta      - udelta(:,2).*sin(theta).*ad_udelta_w(:,2);
ad_udelta_w(:,2)=0;
%12 tl_udelta_w(:,1) = ...
%     + tl_udelta(:,1).*cos(theta) ...
%     - udelta(:,1).*sin(theta).*tl_theta ...
%     + tl_udelta(:,2).*sin(theta) ...
%     + udelta(:,2).*cos(theta).*tl_theta;
ad_udelta(:,1)=ad_udelta(:,1)+ cos(theta)             .*ad_udelta_w(:,1);
ad_theta      =ad_theta      - udelta(:,1).*sin(theta).*ad_udelta_w(:,1);
ad_udelta(:,2)=ad_udelta(:,2)+ sin(theta)             .*ad_udelta_w(:,1);
ad_theta      =ad_theta      + udelta(:,2).*cos(theta).*ad_udelta_w(:,1);
ad_udelta_w(:,1)=0;

% Reniers et al. (2004) model for velocity at top of boundary layer
for i=nx:-1:1
  if(Dr(i)==0)
    %11a2 tl_delta_bl(i)=0;
    ad_delta_bl(i)=0;
    %11a1 tl_udelta(i,:)=[0 0];
    ad_udelta(i,:)=[0 0];
  else
    %11b1 [tl_udelta(i,:),tl_delta_bl(i)] = ...
    %     tl_udelta_reniers2004(tl_ubar(i,:),tl_k(i,:),tl_omega,...
    %                           tl_h(i),tl_Hrms(i),tl_detady(i),...
    %                           tl_tau_wind(i,:),tl_Dr(i),tl_params.fv,tl_params.ks,tl_d50(i),...
    %                           udel_bkgd(i));
    [ad1_ubar,ad1_k,ad1_omega,ad1_h,ad1_Hrms,...
     ad1_detady,ad1_tau_wind,ad1_Dr,ad1_params_fv,ad1_params_ks,ad1_d50] = ...
        ad_udelta_reniers2004(ad_udelta(i,:),ad_delta_bl(i),udel_bkgd(i));

    ad_ubar(i,:)    =ad_ubar(i,:)    +ad1_ubar     ;
    ad_k(i,:)       =ad_k(i,:)       +ad1_k        ;
    ad_omega        =ad_omega        +ad1_omega    ;    
    ad_h(i)         =ad_h(i)         +ad1_h        ;
    ad_Hrms(i)      =ad_Hrms(i)      +ad1_Hrms     ;
    ad_detady(i)    =ad_detady(i)    +ad1_detady   ;
    ad_tau_wind(i,:)=ad_tau_wind(i,:)+ad1_tau_wind ;
    ad_Dr(i)        =ad_Dr(i)        +ad1_Dr       ;
    ad_params.fv    =ad_params.fv    +ad1_params_fv;    
    ad_params.ks    =ad_params.ks    +ad1_params_ks;    
    ad_d50(i)       =ad_d50(i)       +ad1_d50      ;

    ad_udelta(i,:)=0;
    ad_delta_bl(i)=0;
  end
end

% OPTIONAL: Dubarbier et al. suggest a modification to the mean velocity
% prior to calculation of undertow (udelta)
if(params.lambda>0)
  xb=params.lambda*2*pi./kabs;
  ad_xb=zeros(nx,1);  % init
  ad_ur=zeros(nx,2);   % init
  %3 tl_ubar = tl_ur;
  ad_ur = ad_ur + ad_ubar;
  ad_ubar=0;
  for i=1:nx
    ind=find(x(i)-xb(i)<=x&x<=x(i));
    if(length(ind)<=1)
      %2a2 tl_ur(i,2) = tl_ubar0(i,2);
      ad_ubar0(i,2)=ad_ubar0(i,2)+ad_ur(i,2);
      ad_ur(i,2)=0;
      %2a1 tl_ur(i,1) = tl_ubar0(i,1);
      ad_ubar0(i,1)=ad_ubar0(i,1)+ad_ur(i,1);
      ad_ur(i,1)=0;
    else
      xx=xb(i)-x(i)+x(ind);
      xx=xx(:);
      term2=sum(xx);
      ad_xx=zeros(length(ind),1);  % init
      ad_term2=0;  % init
      for j=1:2
        term1=sum(xx.*ubar0(ind,j));
        ad_term1=0;  % init
        %2b7 tl_ur(i,j) = tl_term1/term2 - term1/term2^2*tl_term2;
        ad_term1=ad_term1+ 1/term2      *ad_ur(i,j);
        ad_term2=ad_term2- term1/term2^2*ad_ur(i,j);
        ad_ur(i,j)=0;
        for n=1:length(ind)
          %2b6 tl_term1 = tl_term1 + tl_xx(n)*ubar0(ind(n),j) + xx(n)*tl_ubar0(ind(n),j);
          ad_xx(n)          =ad_xx(n)          + ubar0(ind(n),j)*ad_term1;
          ad_ubar0(ind(n),j)=ad_ubar0(ind(n),j)+ xx(n)          *ad_term1;
        end
        %2b5 tl_term1=0;
        ad_term1=0;
      end
      for n=1:length(ind)
        %2b4 tl_term2 = tl_term2 + tl_xx(n);
        ad_xx(n)=ad_xx(n)+ad_term2;
      end
      %2b3 tl_term2=0;
      ad_term2=0;
      %2b2 tl_xx=tl_xx(:);
      for n=1:length(ind)
        %2b1 tl_xx(n)=tl_xb(i);
        ad_xb(i)=ad_xb(i)+ad_xx(n);
        ad_xx(n)=0;
      end
    end
  end
  % tl_xb = -params.lambda*2*pi./kabs.^2.*tl_kabs ...
  %         + tl_params.lambda*2*pi./kabs;
  ad_kabs=ad_kabs- params.lambda*2*pi./kabs.^2.*ad_xb;
  ad_params.lambda = ad_params.lambda + sum(2*pi./kabs.*ad_xb);
  ad_xb=0;
else
  % tl_ubar = tl_ubar0;
  ad_ubar0 = ad_ubar0 + ad_ubar;
  ad_ubar=0;
end

% settling velocity: use Brown & Lawler.  For vanderA use 0.8*d50
% here, per explanation on page 29
if(strcmp(sedmodel,'vanderA'))
  %10a2 tl_ws = tl_ws_brownLawler(tl_d50_8,.8*d50);
  ad_d50_8 = tl_ws_brownLawler(ad_ws,.8*d50);
  %10a1 tl_d50_8 = .8*tl_d50;
  ad_d50=ad_d50+ .8*ad_d50_8;
  ad_d50_8=0;
else
  %10b1 tl_ws = tl_ws_brownLawler(tl_d50,d50);
  ad_d50 = ad_d50 + tl_ws_brownLawler(ad_ws,d50);
end

% depth averaged mean flow, Nx2 vector, +'ve onshore
%9 tl_ubar0(:,2) = tl_vbar;
ad_vbar=ad_vbar+ad_ubar0(:,2);
ad_ubar0(:,2)=0;
%8 tl_ubar0(:,1) = tl_ubarx;
ad_ubarx=ad_ubarx+ad_ubar0(:,1);
ad_ubar0(:,1)=0;
%7 tl_ubarx = ...
%     - (tl_Ew+2*tl_Er)./(rho*c.*h) ...
%     + (Ew+2*Er)./(rho*c.*h).^2.*rho.*( ...
%         tl_c.*h + tl_h.*c );
ad_Ew=ad_Ew- 1./(rho*c.*h)                   .*ad_ubarx;
ad_Er=ad_Er- 2./(rho*c.*h)                   .*ad_ubarx;
ad_c =ad_c + (Ew+2*Er)./(rho*c.*h).^2.*rho.*h.*ad_ubarx;
ad_h =ad_h + (Ew+2*Er)./(rho*c.*h).^2.*rho.*c.*ad_ubarx;
ad_ubarx=0;

% convert from (kabs,theta) to vector wavenumber, and calculate c for
% convenience
%6 tl_c = ...
%     + tl_omega./kabs ...
%     - omega./kabs.^2.*tl_kabs;
ad_omega=ad_omega+ sum(1./kabs       .*ad_c);
ad_kabs =ad_kabs - omega./kabs.^2.*ad_c;
ad_c=0;
%5 tl_k(:,2) = ...
%     + tl_kabs(:).*sin(theta(:)) ...
%     + kabs(:).*cos(theta(:)).*tl_theta;
ad_kabs(:)=ad_kabs(:)+ sin(theta(:))         .*ad_k(:,2);
ad_theta  =ad_theta  + kabs(:).*cos(theta(:)).*ad_k(:,2);
ad_k(:,2)=0;
%4 tl_k(:,1) = ...
%     + tl_kabs(:).*cos(theta(:)) ...
%     - kabs(:).*sin(theta(:)).*tl_theta;
ad_kabs(:)=ad_kabs(:)+ cos(theta(:))         .*ad_k(:,1);
ad_theta  =ad_theta  - kabs(:).*sin(theta(:)).*ad_k(:,1);
ad_k(:,1)=0;

% special case dt=0: in this case we are done since there is no morphology
% update to compute
if(dt==0)
  %3a2 tl_hp=tl_h;
  ad_h = ad_h + ad_hp;
  ad_hp=0;
  %3a1 tl_Qx=zeros(nx,1);
  ad_Qx=zeros(nx,1);
else

% wave shape parameters.  Note Uwave_ruessink2012 specifies Hmo as input
% tl_Sw=tl_Sw0+tl_dSw;
ad_Sw0=ad_Sw0+ad_Sw;
ad_dSw=ad_dSw+ad_Sw;
ad_Sw=zeros(nx,1);
% tl_Aw=tl_Aw0+tl_dAw;
ad_Aw0=ad_Aw0+ad_Aw;
ad_dAw=ad_dAw+ad_Aw;
ad_Aw=zeros(nx,1);
% [tl_Aw0,tl_Sw0,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_kabs,tl_omega,tl_h,uwave_bkgd);
[ad1_Hmo,ad1_kabs,ad1_omega,ad1_h]=ad_Uwave_ruessink2012_params(ad_Aw0,ad_Sw0,ad_Uw,uwave_bkgd);
ad_Hmo  =ad_Hmo  +ad1_Hmo  ;
ad_kabs =ad_kabs +ad1_kabs ;
ad_omega=ad_omega+ad1_omega;
ad_h    =ad_h    +ad1_h    ;
% tl_Hmo=1.4*tl_Hrms;
ad_Hrms=ad_Hrms+1.4*ad_Hmo;
ad_Hmo=0;

% 1DH wave and longshore current balance
% [tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
%     tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,tl_detady,tl_dgamma,tl_beta0,hydro_bkgd);
[ad1_h,ad1_H0,ad1_theta0,ad1_omega,ad1_ka_drag,ad1_tau_wind,ad1_detady,ad1_dgamma,ad1_beta0] = ...
    ad_hydro_ruessink2001(ad_Hrms,ad_theta,ad_vbar,ad_kabs,ad_Ew,ad_Er,ad_Dr,hydro_bkgd);
ad_h       =ad_h       +ad1_h       ;
ad_H0      =ad_H0      +ad1_H0      ;
ad_theta0  =ad_theta0  +ad1_theta0  ;
ad_omega   =ad_omega   +ad1_omega   ;
ad_ka_drag =ad_ka_drag +ad1_ka_drag ;
ad_tau_wind=ad_tau_wind+ad1_tau_wind;
ad_detady  =ad_detady  +ad1_detady  ;
ad_dgamma  =ad_dgamma  +ad1_dgamma  ;
ad_beta0   =ad_beta0   +ad1_beta0   ;

%1 tl_h(imask)=0;  % min depth constraint
ad_h(imask)=0;

end  % catch for special case dt==0

end  % main function for single time step
