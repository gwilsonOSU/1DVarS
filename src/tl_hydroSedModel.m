function [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hpout] = ...
    tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                     tl_detady,tl_dgamma,tl_dAw,tl_dSw,...
                     tl_d50,tl_d90,tl_params,bkgd)%,outvar)

% sub-stepping loop
tl_hp(:,1) = tl_h;  % init t=0
for n=1:bkgd(1).nsubsteps
  % disp(['hydroSedModel substep ' num2str(n) ' of ' num2str(bkgd(1).nsubsteps)])
  [tl_Hrms(:,n),tl_vbar(:,n),tl_theta(:,n),tl_kabs(:,n),...
   tl_Qx(:,n),tl_hp(:,n+1)] = ...
      tl_hydroSedModel_main(tl_hp(:,n),tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,tl_detady,...
                         tl_dgamma,tl_dAw,tl_dSw,tl_d50,tl_d90,tl_params,bkgd(n));%,outvar);
end

% drop initial condition from hpout
for n=1:bkgd(1).nsubsteps
  tl_hpout(:,n)=tl_hp(:,n+1);
end

end  % end wrapper function (for sub-stepping loop logic)

% begin main function, for single time step
function [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hp] = ...
    tl_hydroSedModel_main(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,tl_detady,tl_dgamma,tl_dAw,tl_dSw,...
                  tl_d50,tl_d90,tl_params,bkgd)%,outvar)

physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

nx=length(x);

%--------------------------------------
% begin TL
%--------------------------------------

tl_h(imask)=0;  % min depth constraint

% 1DH wave and longshore current balance
[tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
    tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,tl_detady,tl_dgamma,tl_beta0,hydro_bkgd);

% wave shape parameters.  Note Uwave_ruessink2012 specifies Hmo as input
tl_Hmo=1.4*tl_Hrms;
[tl_Aw0,tl_Sw0,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_kabs,tl_omega,tl_h,uwave_bkgd);
tl_Aw=tl_Aw0+tl_dAw;
tl_Sw=tl_Sw0+tl_dSw;

% special case dt=0: in this case we are done since there is no morphology
% update to compute
if(dt==0)
  tl_Qx=zeros(nx,1);
  tl_hp=tl_h;
else

% convert from (kabs,theta) to vector wavenumber, and calculate c for
% convenience
tl_k(:,1) = ...
    + tl_kabs(:).*cos(theta(:)) ...
    - kabs(:).*sin(theta(:)).*tl_theta;
tl_k(:,2) = ...
    + tl_kabs(:).*sin(theta(:)) ...
    + kabs(:).*cos(theta(:)).*tl_theta;
tl_c = ...
    + tl_omega./kabs ...
    - omega./kabs.^2.*tl_kabs;

% depth averaged mean flow, Nx2 vector, +'ve onshore
tl_ubarx = ...
    - (tl_Ew+2*tl_Er)./(rho*c.*h) ...
    + (Ew+2*Er)./(rho*c.*h).^2.*rho.*( ...
        tl_c.*h + tl_h.*c );
tl_ubar0(:,1) = tl_ubarx;
tl_ubar0(:,2) = tl_vbar;

% settling velocity: use Brown & Lawler.  TODO, for vanderA use 0.8*d50
% here, per explanation on page 29
if(strcmp(sedmodel,'vanderA'))
  tl_d50_8 = .8*tl_d50;
  tl_ws = tl_ws_brownLawler(tl_d50_8,.8*d50);
else
  tl_ws = tl_ws_brownLawler(tl_d50,d50);
end

% OPTIONAL: Dubarbier et al. suggest a modification to the mean velocity
% prior to calculation of undertow (udelta)
if(params.lambda>0)
  xb=params.lambda*2*pi./kabs;
  tl_xb = -params.lambda*2*pi./kabs.^2.*tl_kabs ...
          + tl_params.lambda*2*pi./kabs;
  for i=1:nx
    ind=find(x(i)-xb(i)<=x&x<=x(i));
    if(length(ind)<=1)
      tl_ur(i,1) = tl_ubar0(i,1);
      tl_ur(i,2) = tl_ubar0(i,2);
    else
      xx=xb(i)-x(i)+x(ind);
      xx=xx(:);
      for n=1:length(ind)
        tl_xx(n)=tl_xb(i);
      end
      tl_xx=tl_xx(:);
      term2=sum(xx);
      tl_term2=0;
      for n=1:length(ind)
        tl_term2 = tl_term2 + tl_xx(n);
      end
      for j=1:2
        term1=sum(xx.*ubar0(ind,j));
        tl_term1=0;
        for n=1:length(ind)
          tl_term1 = tl_term1 + tl_xx(n)*ubar0(ind(n),j) + xx(n)*tl_ubar0(ind(n),j);
        end
        tl_ur(i,j) = tl_term1/term2 - term1/term2^2*tl_term2;
      end
    end
  end
  tl_ubar = tl_ur;
else
  tl_ubar = tl_ubar0;
end

% Reniers et al. (2004) model for velocity at top of boundary layer
tl_delta_bl=zeros(nx,1);  % init
tl_udelta=tl_ubar;   % init
for i=1:nx
  if(Dr(i)>0)
    [tl_udelta(i,:),tl_delta_bl(i)] = ...
        tl_udelta_reniers2004(tl_ubar(i,:),tl_k(i,:),tl_omega,...
                              tl_h(i),tl_Hrms(i),tl_detady(i),...
                              tl_tau_wind(i,:),tl_Dr(i),tl_params.fv,tl_params.ks,tl_d50(i),...
                              udel_bkgd(i));
  end
end

% rotate udelta into wave direction, as assumed by sed transport equations
% udelta_w(:,1) = +udelta(:,1).*cos(theta) + udelta(:,2).*sin(theta);
tl_udelta_w(:,1) = ...
    + tl_udelta(:,1).*cos(theta) ...
    - udelta(:,1).*sin(theta).*tl_theta ...
    + tl_udelta(:,2).*sin(theta) ...
    + udelta(:,2).*cos(theta).*tl_theta;
% udelta_w(:,2) = -udelta(:,1).*sin(theta) + udelta(:,2).*cos(theta);
tl_udelta_w(:,2) = ...
    - tl_udelta(:,1).*sin(theta) ...
    - udelta(:,1).*cos(theta).*tl_theta ...
    + tl_udelta(:,2).*cos(theta) ...
    - udelta(:,2).*sin(theta).*tl_theta;

% run the requested model for sediment flux (m2/s)
tl_tanbeta = tl_calcTanbeta(tl_h,x)';
if(strcmp(sedmodel,'dubarbier'))  % Dubarbier et al. (2015)
  tl_Q0 =tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,...
                            tl_params.Cw,tl_params.Cc,tl_params.Cf,tl_params.Ka,...
                            bkgd_qtrans);
elseif(strcmp(sedmodel,'soulsbyVanRijn'))  % Soulsby & van Rijn
  tl_Q0 = tl_qtrans_soulsbyVanRijn(tl_d50,tl_d90,tl_h,tl_tanbeta,...
                                  tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubar,...
                                  tl_Dr,tl_param,bkgd_qtrans);
elseif(strcmp(sedmodel,'vanderA'))  % van Der A et al. (2013)
  tl_Q0 = tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,...
                           tl_udelta,tl_delta_bl,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_params,bkgd_qtrans);
end

% mitigate transport discontinuity at the shoreline
% Q1=Q0;
tl_Q1 = tl_Q0;
% Q1(imax:end)=0;
tl_Q1(imax:end)=0;

% filtering is applied in NL model, neglected in TL
% Q = Q1;
tl_Q = tl_Q1;

% rotate output from wave-following coords to cartesian
tl_Qx = ...
    + tl_Q.*cos(theta) ...
    - Q.*sin(theta).*tl_theta;

% bathymetry update: dhdt = -dzdt = dQdx.  This is the Exner equation,
% e.g. see Dubarbier et al. (2015) eqn. (16), and note Q is the volumetric
% transport rate (m2/s) including the bed porosity
if(doMarieu)  % use "stable" Marieu formulation for dh/dt
  dx=abs(x(2)-x(1));
  [tl_hp1,tl_qp1]=tl_dhdt_marieu2007(tl_Qx ,tl_h  ,bkgd_marieu_step1);
  [tl_hp ,tl_qp ]=tl_dhdt_marieu2007(tl_qp1,tl_hp1,bkgd_marieu_step2);
  tl_dh = tl_hp - tl_h;
  tl_dQdx=zeros(nx,1);
else
  tl_dQdx = tl_ddx_upwind(tl_Qx,x,Qx,horig);
  % tl_dQdx = tl_ddx_centered(tl_Qx,x,Qx);
  tl_dh = tl_dQdx*dt;
  tl_qp = zeros(nx,1);
end
tl_dh = tl_dh.*wgt;   % apply damping near shore
tl_dh(isnan(dh))=0;
tl_hp = tl_h + tl_dh;

end  % catch for special case dt==0

% % TEST-CODE: override output variable
% if((length(outvar)>=6 & strcmp(outvar(1:6),'udelta')) | strcmp(outvar,'ubar'))
%   eval(['tl_hp = tl_' outvar '(:,2);']);
%   eval(['tl_' outvar '(:,1)=0;']);
% else
%   eval(['tl_hp = tl_' outvar ';']);
% end
% tl_vbar =zeros(nx,1);
% tl_theta=zeros(nx,1);
% tl_kabs =zeros(nx,1);
% tl_Hrms =zeros(nx,1);
% tl_Qx   =zeros(nx,1);

% TEST: tweak output for TL testing
% eval(['tl_Qx=tl_' outvar ';']);

end  % main function for single time step
