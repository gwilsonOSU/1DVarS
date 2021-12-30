%
% Full-run adjoint sensitivity. The final-time adjoint bathymetry is
% initialized with the identity matrix and propagates backwards through
% time. Adjoints for time-independent parameters are stored for each time
% step, and should be summed to get the overall adjoint sensitivity for
% final bathymetry.
%
% Examples:
%
% ex1) ad_theta0(i,n) represents the sensitivity of final bathymetry at
%       gridpoint i to unit perturbations in theta0 at time n.
%
% ex2) ad_h(i,j,n) represents the sensitivity of final bathymetry at gridpoint
%      i to unit perturbations in bathymetry at gridpoint j at time n.
%
% ex2b) The matrix ad_h(:,:,1) shows how initial bathymetry affects final
%        bathymetry.
%
% TODO: Instead of final-bathymetry, this code could easily be adapted to
% compute adjoint sensitivity for sediment flux (Qx), or other variables.
%
addpath util
addpath(genpath('../../src'))
clear all

% USER-INPUT: Choose a case
duck94Case='b';

parpoolN=2;

%-------------------------------------------------
% End of user input.  No edits needed beyond this point.
%-------------------------------------------------

% Load the model grid and grid-referenced observational data for this case.
% Cached as mat-files to save time.
obsdatafn=['obsdataCache/obsdata_case' duck94Case '.mat'];
if(~isempty(dir(obsdatafn)))
  disp(['loading pre-cached obsdata function: ' obsdatafn])
  load(obsdatafn);
else
  disp(['loading obsdata'])
  [hydroobs,bathyobs,grid,waves8m,windEOP]=prepObsData(dnum,bathyfn,duck94Case);
  disp(['caching obsdata for next time: ' obsdatafn])
  save(obsdatafn,'hydroobs','bathyobs','grid','waves8m','windEOP');
end

% start a parallel pool
currentPool=gcp('nocreate');
if(isempty(currentPool) | currentPool.NumWorkers ~= parpoolN)
  if(~isempty(currentPool))
    delete(gcp('nocreate'));
  end
  parpool('local',parpoolN);
end

% use saved outputs from duck94TwoPhaseInversion.m as a background state for
% the analysis.  Do a full model run to regenerate the time-dependent
% background.
load(['case_' duck94Case '_outputs/assimIter1/output.mat']);
modelinput=initModelInputs(duck94Case,grid,bkgd_1.sedmodel);
modelinput.params=bkgd_1.params;
hydroobs=hydroobs(1:bathyobs(end).obsn);  % end at bathyobs time
bkgd=hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs);
nx=grid.nx;
nt=length(bkgd);

disp('continue manually')
return;

%-------------------------------------------------------
% V1, test code: Full-run adjoint for h
%-------------------------------------------------------

% init full-run adjoint variables
ad_full_tau_windx=zeros(nx,nx,nt);
ad_full_tau_windy=zeros(nx,nx,nt);
ad_full_h        =zeros(nx,nx,nt);
ad_full_H0       =zeros(nx,nt);
ad_full_theta0   =zeros(nx,nt);
ad_full_omega    =zeros(nx,nt);
ad_full_ka_drag  =zeros(nx,1);
ad_full_detady   =zeros(nx,nx,nt);
ad_full_dgamma   =zeros(nx,nx,nt);
ad_full_dAw      =zeros(nx,nx,nt);
ad_full_dSw      =zeros(nx,nx,nt);
ad_full_d50      =zeros(nx,nx,nt);
ad_full_d90      =zeros(nx,nx,nt);
ad_full_params.fv=zeros(nx,1);
ad_full_params.ks=zeros(nx,1);
if(strcmp(bkgd(1).sedmodel,'vanderA'))
  ad_full_params.n    =zeros(nx,1);
  ad_full_params.m    =zeros(nx,1);
  ad_full_params.xi   =zeros(nx,1);
  ad_full_params.alpha=zeros(nx,1);
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  ad_full_params.Cw=zeros(nx,1);
  ad_full_params.Cc=zeros(nx,1);
  ad_full_params.Cf=zeros(nx,1);
  ad_full_params.Ka=zeros(nx,1);
end

% calculate full-run bathymetry adjoints at every gridpoint
for i=1:nx  % TODO, use parfor, will require disk caching though

  % initialize comb at gridpoint i and final time step
  ad_Hrms =zeros(nx,1);
  ad_theta=zeros(nx,1);
  ad_vbar =zeros(nx,1);
  ad_kabs =zeros(nx,1);
  ad_Qx   =zeros(nx,1);
  ad_h   =zeros(nx,nt+1);
  ad_h(i,nt+1)=1;  % comb

  % initialize adjoint outputs
  ad_params=paramsHandler(0,bkgd(1).sedmodel,0,0,0,0,0,0);  % init ad_params struct to zero
  ad_ka_drag=0;
  ad_d50=zeros(nx,1);
  ad_d90=zeros(nx,1);
  ad_H0=zeros(nt,1);
  ad_theta0  =zeros(nt,1);
  ad_omega   =zeros(nt,1);
  ad_dgamma  =zeros(nx,nt);
  ad_dAw     =zeros(nx,nt);
  ad_dSw     =zeros(nx,nt);
  ad_tau_wind=zeros(nx,2,nt);
  ad_detady  =zeros(nx,nt);

  % propagate adjoint backwards from time nt to 1
  for n2=nt:-1:1
    [ad_h(:,n2),ad_H0(n2),ad_theta0(n2),ad_omega(n2),ad1_ka_drag,ad_tau_wind(:,:,n2),...
     ad_detady(:,n2),ad_dgamma(:,n2),ad_dAw(:,n2),ad_dSw(:,n2),...
     ad1_d50,ad1_d90,ad1_params] = ...
        ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_h(:,n2+1),bkgd(n2));
    ad_ka_drag     =ad_ka_drag     +ad1_ka_drag     ;
    ad_d50         =ad_d50         +ad1_d50         ;
    ad_d90         =ad_d90         +ad1_d90         ;
    ad_params.fv   =ad_params.fv   +ad1_params.fv   ;
    ad_params.ks   =ad_params.ks   +ad1_params.ks   ;
    if(strcmp(bkgd(1).sedmodel,'vanderA'))
      ad_params.n    =ad_params.n    +ad1_params.n    ;
      ad_params.m    =ad_params.m    +ad1_params.m    ;
      ad_params.xi   =ad_params.xi   +ad1_params.xi   ;
      ad_params.alpha=ad_params.alpha+ad1_params.alpha;
    elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
      ad_params.Cw = ad_params.Cw + ad1_params.Cw;
      ad_params.Cc = ad_params.Cc + ad1_params.Cc;
      ad_params.Cf = ad_params.Cf + ad1_params.Cf;
      ad_params.Ka = ad_params.Ka + ad1_params.Ka;
    end
  end

  % save adjoint info for this grid point
  ad_full_tau_windx(i,:,:)=squeeze(ad_tau_wind(:,1,:));
  ad_full_tau_windy(i,:,:)=squeeze(ad_tau_wind(:,2,:));
  ad_full_h        (i,:,:)=ad_h        ;
  ad_full_H0       (i,:)  =ad_H0       ;
  ad_full_theta0   (i,:)  =ad_theta0   ;
  ad_full_omega    (i,:)  =ad_omega    ;
  ad_full_ka_drag  (i)    =ad_ka_drag  ;
  ad_full_detady   (i,:,:)=ad_detady   ;
  ad_full_dgamma   (i,:,:)=ad_dgamma   ;
  ad_full_dAw      (i,:,:)=ad_dAw      ;
  ad_full_dSw      (i,:,:)=ad_dSw      ;
  ad_full_d50      (i,:,:)=ad_d50      ;
  ad_full_d90      (i,:,:)=ad_d90      ;
  ad_full_params.fv(i)    =ad_params.fv;
  ad_full_params.ks(i)    =ad_params.ks;
  if(strcmp(bkgd(1).sedmodel,'vanderA'))
    ad_full_params.n    (i)=ad_full_params.n    ;
    ad_full_params.m    (i)=ad_full_params.m    ;
    ad_full_params.xi   (i)=ad_full_params.xi   ;
    ad_full_params.alpha(i)=ad_full_params.alpha;
  elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
    ad_full_params.Cw(i)=ad_params.Cw;
    ad_full_params.Cc(i)=ad_params.Cc;
    ad_full_params.Cf(i)=ad_params.Cf;
    ad_full_params.Ka(i)=ad_params.Ka;
  end

end  % loop for adjoint at each gridpoint

% rename the "full" adjoint variables for simplicity
ad_tau_windx_wind=ad_full_tau_windx_wind;
ad_tau_windy_wind=ad_full_tau_windy_wind;
ad_h             =ad_full_h             ;
ad_H0            =ad_full_H0            ;
ad_theta0        =ad_full_theta0        ;
ad_omega         =ad_full_omega         ;
ad_ka_drag       =ad_full_ka_drag       ;
ad_detady        =ad_full_detady        ;
ad_dgamma        =ad_full_dgamma        ;
ad_dAw           =ad_full_dAw           ;
ad_dSw           =ad_full_dSw           ;
ad_d50           =ad_full_d50           ;
ad_d90           =ad_full_d90           ;
ad_params.fv     =ad_full_params.fv     ;
ad_params.ks     =ad_full_params.ks     ;
if(strcmp(bkgd(1).sedmodel,'vanderA'))
  ad_params.n    =ad_full_params.n    ;
  ad_params.m    =ad_full_params.m    ;
  ad_params.xi   =ad_full_params.xi   ;
  ad_params.alpha=ad_full_params.alpha;
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  ad_params.Cw=ad_full_params.Cw;
  ad_params.Cc=ad_full_params.Cc;
  ad_params.Cf=ad_full_params.Cf;
  ad_params.Ka=ad_full_params.Ka;
end
clear ad_full_*

%-------------------------------------------------------
% V2, test code: Single-time adjoint for Qx
%
% example: ad_d50(j,i) is the sensitivity of Qx(i) to perturbations in
% d50(j).  Can sum over the 1st dimension to get overall sensitivity to d50
% at all gridpoints.
%-------------------------------------------------------

% case-dependent parameters. Note 'n' is the time step to be analyzed, need
% to pick a representative one.
if(duck94Case=='b')
  n=200;  % 200 is a high-Qx time for case-b
  qscale=2e-4;
  qsign=+1;
elseif(duck94Case=='c')
  n=250;  % 250 is a high-Qx time for case-c
  qscale=1e-3;
  qsign=+1;
else
  disp('choose case manually')
  return;
end

% init AD outputs
ad_h=nan(nx,nx);
ad_H0=nan(nx,1);
ad_theta0=nan(nx,1);
ad_omega=nan(nx,1);
ad_ka_drag=nan(nx,1);
ad_tau_wind=nan(nx,2,nx);
ad_detady=nan(nx,nx);
ad_dgamma=nan(nx,nx);
ad_dAw=nan(nx,nx);
ad_dSw=nan(nx,nx);
ad_d50=nan(nx,nx);
ad_d90=nan(nx,nx);

for i=1:nx  % TODO, use parfor?

  % initialize comb at gridpoint i
  ad_Hrms =zeros(nx,1);
  ad_theta=zeros(nx,1);
  ad_vbar =zeros(nx,1);
  ad_kabs =zeros(nx,1);
  ad_Qx   =zeros(nx,1);
  ad_hp   =zeros(nx,1);
  ad_Qx(i)=1;  % comb

  % run adjoint for time n
  [ad_h(:,i),ad_H0(i),ad_theta0(i),ad_omega(i),ad_ka_drag(i),ad_tau_wind(:,:,i),...
   ad_detady(:,i),ad_dgamma(:,i),ad_dAw(:,i),ad_dSw(:,i),ad_d50(:,i),ad_d90(:,i),ad_params(i)] = ...
      ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hp,bkgd(n));

end  % loop for adjoint at each gridpoint

% normalize to Qx-units by mulitiplying by the bkgd value of each variable,
% and apply qsign
for i=1:nx
  ad_h_norm     (:,i)   =ad_h     (:,i)    .*bkgd(n).h               *qsign; 
  ad_detady_norm(:,i)   =ad_detady(:,i)    .*bkgd(n).detady          *qsign;
  ad_gamma_norm (:,i)   =ad_dgamma(:,i)    .*bkgd(n).hydro_bkgd.gamma*qsign;
  ad_Aw_norm    (:,i)   =ad_dAw   (:,i)    .*bkgd(n).Aw              *qsign;
  ad_Sw_norm    (:,i)   =ad_dSw   (:,i)    .*bkgd(n).Sw              *qsign;
  ad_d50_norm   (:,i)   =ad_d50   (:,i)    .*bkgd(n).d50             *qsign;
  ad_d90_norm   (:,i)   =ad_d90   (:,i)    .*bkgd(n).d90             *qsign;
  ad_tau_windx_norm(:,i)=ad_tau_wind(:,1,i).*bkgd(n).tau_wind(:,1)   *qsign;
  ad_tau_windy_norm(:,i)=ad_tau_wind(:,2,i).*bkgd(n).tau_wind(:,2)   *qsign;
end
ad_fv_norm   =[ad_params.fv   ]*bkgd(n).params.fv   *qsign;
ad_ks_norm   =[ad_params.ks   ]*bkgd(n).params.ks   *qsign;
ad_n_norm    =[ad_params.n    ]*bkgd(n).params.n    *qsign;
ad_m_norm    =[ad_params.m    ]*bkgd(n).params.m    *qsign;
ad_xi_norm   =[ad_params.xi   ]*bkgd(n).params.xi   *qsign;
ad_alpha_norm=[ad_params.alpha]*bkgd(n).params.alpha*qsign;
ad_H0_norm     =ad_H0     *bkgd(n).H0     *qsign;
ad_theta0_norm =ad_theta0 *bkgd(n).theta0 *qsign;
ad_omega_norm  =ad_omega  *bkgd(n).omega  *qsign;
ad_ka_drag_norm=ad_ka_drag*bkgd(n).ka_drag*qsign;

% show results for vector variables.  Skip tau_wind and detady, it makes
% sense that they have negligible sensitivity
figure(1),clf
subplot(2,2,1); pcolor(grid.xFRF,grid.xFRF,ad_h_norm        ),sf,ylabel('\Delta h @ x')
subplot(2,2,2); pcolor(grid.xFRF,grid.xFRF,ad_gamma_norm    ),sf,ylabel('\Delta \gamma @ x')
subplot(2,2,3); pcolor(grid.xFRF,grid.xFRF,ad_Aw_norm       ),sf,ylabel('\Delta Aw @ x')
subplot(2,2,4); pcolor(grid.xFRF,grid.xFRF,ad_Sw_norm       ),sf,ylabel('\Delta Sw @ x')
for i=1:4
  subplot(2,2,i)
  caxis([-1 1]*qscale)
  axis([100 300 100 300])
  xlabel('\Delta Q @ x')
end

% Overlay all sensitivities, grouped by variable type
figure(2),clf
subplot(2,2,1)  % sed params
hold on
lstr={};
title('Sediment Transport Parameters')
plot(grid.xFRF,ad_n_norm      ),lstr{end+1}='\Delta n       ';
plot(grid.xFRF,ad_m_norm      ),lstr{end+1}='\Delta m       ';
plot(grid.xFRF,ad_xi_norm     ),lstr{end+1}='\Delta \xi   ';
plot(grid.xFRF,ad_alpha_norm  ),lstr{end+1}='\Delta \alpha';
legend(lstr)
subplot(2,2,2)  % hydro BCs
hold on
lstr={};
title('Wave Model BCs')
plot(grid.xFRF,-ad_H0_norm    ),lstr{end+1}='-\Delta H_0     ';
plot(grid.xFRF,ad_omega_norm  ),lstr{end+1}='\Delta \omega  ';
plot(grid.xFRF,ad_theta0_norm ),lstr{end+1}='\Delta \theta_0'; 
plot(grid.xFRF,max(ad_h_norm )),lstr{end+1}='\Delta h     '; 
legend(lstr)
subplot(2,2,3)  % hydro interior errors (+bathymetry)
hold on
lstr={};
title('Wave Model Interior Variables')
plot(grid.xFRF,max(ad_gamma_norm)),lstr{end+1}='\Delta \gamma'; 
plot(grid.xFRF,max(ad_Aw_norm   )),lstr{end+1}='\Delta Aw    '; 
plot(grid.xFRF,max(ad_Sw_norm   )),lstr{end+1}='\Delta Sw    '; 
legend(lstr)
subplot(2,2,4)  % undertow params
hold on
lstr={};
title('Undertow Model Parameters')
plot(grid.xFRF,ad_fv_norm),lstr{end+1}='\Delta f_v     ';
plot(grid.xFRF,ad_ks_norm),lstr{end+1}='\Delta k_s   ';
legend(lstr)
for i=1:4
  subplot(2,2,i)
  box on
  ylim([-1 1]*qscale)
  xlim([150 400])
  ylabel('\Delta Q [m^2/s]')
  xlabel('x [m]')
end
