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
clearvars -except duck94Case

% USER-INPUT: Choose a case
% duck94Case='b';

cacheDir='/tmp/bathyAssimCache';

doprint=1;

% pick a representative time step 'targetn' for adjoint analysis, and define
% some plotting parameters
if(duck94Case=='b')
  targetn=200;  % 200 is a high-Qx time for case-b
  qscale=2e-4;
  qsign=+1;
elseif(duck94Case=='c')
  targetn=250;  % 250 is a high-Qx time for case-c
  qscale=1e-3;
  qsign=+1;
else
  disp('choose case manually')
  return;
end

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

% % use saved outputs from duck94TwoPhaseInversion.m as a background state for
% % the analysis.  Do a full model run to regenerate the time-dependent
% % background.
% load(['case_' duck94Case '_outputs/assimIter1/output.mat']);
% modelinput=initModelInputs(duck94Case,grid,bkgd_1.sedmodel);
% modelinput.params=bkgd_1.params;
% hydroobs=hydroobs(1:bathyobs(end).obsn);  % end at bathyobs time

% use default model inputs
modelinput=initModelInputs(duck94Case,grid,sedmodel);

% run forward model (with hydro assimilation)
hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs,cacheDir);

disp('continue manually')
return;

%-------------------------------------------------------
% V2, test code: Single-time adjoint for Qx
%
% example: ad_d50(j,i) is the sensitivity of Qx(i) to perturbations in
% d50(j).  Can sum over the 1st dimension to get overall sensitivity to d50
% at all gridpoints.
%-------------------------------------------------------

n=targetn;
nx=grid.nx;
nt=length(hydroobs);

% init AD outputs
ad_h=nan(nx,nx);
ad_H0=nan(nx,1);
ad_theta0=nan(nx,1);
ad_omega=nan(nx,1);
ad_ka_drag=nan(nx,1);
ad_beta0  =nan(nx,1);
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
  bkgdn=load([cacheDir '/bkgd' num2str(n) '.mat']);
  [ad_h(:,i),ad_H0(i),ad_theta0(i),ad_omega(i),ad_ka_drag(i),ad_beta0(i),ad_tau_wind(:,:,i),...
   ad_detady(:,i),ad_dgamma(:,i),ad_dAw(:,i),ad_dSw(:,i),ad_d50(:,i),ad_d90(:,i),ad_params(i)] = ...
      ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_hp,bkgdn);

end  % loop for adjoint at each gridpoint

% normalize to Qx-units by mulitiplying by the bkgd value of each variable,
% and apply qsign
for i=1:nx
  ad_h_norm     (:,i)   =ad_h     (:,i)    .*bkgdn.h               *qsign; 
  ad_detady_norm(:,i)   =ad_detady(:,i)    .*bkgdn.detady          *qsign;
  ad_gamma_norm (:,i)   =ad_dgamma(:,i)    .*bkgdn.hydro_bkgd.gamma*qsign;
  ad_Aw_norm    (:,i)   =ad_dAw   (:,i)    .*bkgdn.Aw              *qsign;
  ad_Sw_norm    (:,i)   =ad_dSw   (:,i)    .*bkgdn.Sw              *qsign;
  ad_d50_norm   (:,i)   =ad_d50   (:,i)    .*bkgdn.d50             *qsign;
  ad_d90_norm   (:,i)   =ad_d90   (:,i)    .*bkgdn.d90             *qsign;
  ad_tau_windx_norm(:,i)=ad_tau_wind(:,1,i).*bkgdn.tau_wind(:,1)   *qsign;
  ad_tau_windy_norm(:,i)=ad_tau_wind(:,2,i).*bkgdn.tau_wind(:,2)   *qsign;
end
ad_fv_norm   =[ad_params.fv   ]*bkgdn.params.fv   *qsign;
ad_ks_norm   =[ad_params.ks   ]*bkgdn.params.ks   *qsign;
ad_lambda_norm=[ad_params.lambda]*bkgdn.params.lambda*qsign;
ad_n_norm    =[ad_params.n    ]*bkgdn.params.n    *qsign;
ad_m_norm    =[ad_params.m    ]*bkgdn.params.m    *qsign;
ad_xi_norm   =[ad_params.xi   ]*bkgdn.params.xi   *qsign;
ad_alpha_norm=[ad_params.alpha]*bkgdn.params.alpha*qsign;
ad_Cc_norm=[ad_params.Cc]*bkgdn.params.Cc*qsign;
ad_Cf_norm=[ad_params.Cf]*bkgdn.params.Cf*qsign;
ad_H0_norm     =ad_H0     *bkgdn.H0     *qsign;
ad_theta0_norm =ad_theta0 *bkgdn.theta0 *qsign;
ad_omega_norm  =ad_omega  *bkgdn.omega  *qsign;
ad_ka_drag_norm=ad_ka_drag*bkgdn.ka_drag*qsign;
ad_beta0_norm  =ad_beta0  *bkgdn.beta0  *qsign;

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
if(doprint)
  print -dpng duck94AdjointSensitivityOneTime_Qtermsvec.png
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
plot(grid.xFRF,ad_Cc_norm),lstr{end+1}='\Delta C_c';
plot(grid.xFRF,ad_Cf_norm),lstr{end+1}='\Delta C_f';
legend(lstr)
subplot(2,2,2)  % hydro BCs
hold on
lstr={};
title('Wave Model BCs')
plot(grid.xFRF,-ad_H0_norm    ),lstr{end+1}='-\Delta H_0     ';
plot(grid.xFRF,ad_omega_norm  ),lstr{end+1}='\Delta \omega  ';
plot(grid.xFRF,ad_theta0_norm ),lstr{end+1}='\Delta \theta_0'; 
plot(grid.xFRF,ad_beta0_norm  ),lstr{end+1}='\Delta \beta_0'; 
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
plot(grid.xFRF,ad_lambda_norm),lstr{end+1}='\Delta \lambda';
legend(lstr)
for i=1:4
  subplot(2,2,i)
  box on
  ylim([-1 1]*qscale)
  xlim([150 400])
  ylabel('\Delta Q [m^2/s]')
  xlabel('x [m]')
end
if(doprint)
  print -dpng duck94AdjointSensitivityOneTime_Qterms.png
end
