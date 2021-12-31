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
% NOTE: This function will use the local directory /tmp/bathyAssimCache as a
% disk cache to work around memory limitation.  If the directory doesn't
% exist, it will be created, and the cached data will be deleted when
% finished..  You must have ~5GB free in /tmp, for typical runs with ~500
% time steps.
%
addpath util
addpath(genpath('../../src'))
clearvars -except duck94Case parpoolN

% USER-INPUT: Choose a case
% duck94Case='b';

parpoolN=12;

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
% Full-run adjoint for h
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

% start a parallel pool
currentPool=gcp('nocreate');
if(isempty(currentPool) | currentPool.NumWorkers ~= parpoolN)
  if(~isempty(currentPool))
    delete(gcp('nocreate'));
  end
  parpool('local',parpoolN);
end

% use disk caching of background struct to work around shared-memory
% limitation of matlab's parallel toolbox
disp('Caching bkgd struct to disk')
tmpdir='/tmp/bathyAssimCache';
if(isempty(dir(tmpdir)))
  mkdir(tmpdir);
end
nt=length(bkgd);
for n=1:nt
  if(floor(n/nt*10)>floor((n-1)/nt*10))
    disp(['   ' num2str(floor(n/nt*100)) '%'])
  end
  this=bkgd(n);
  save([tmpdir '/bkgd' num2str(n) '.mat'],'-struct','this');
end

% calculate full-run bathymetry adjoints at every gridpoint
parfor i=1:nx

  % initialize comb at gridpoint i and final time step
  ad_Hrms =zeros(nx,1);
  ad_theta=zeros(nx,1);
  ad_vbar =zeros(nx,1);
  ad_kabs =zeros(nx,1);
  ad_Qx   =zeros(nx,1);
  ad_h   =zeros(nx,nt+1);
  ad_h(i,nt+1)=1;  % comb

  % initialize adjoint outputs
  bkgd1=load([tmpdir '/bkgd1.mat']);
  ad_params=paramsHandler(0,bkgd1.sedmodel,0,0,0,0,0,0);  % init ad_params struct to zero
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
    bkgdn2=load([tmpdir '/bkgd' num2str(n2) '.mat']);
    [ad_h(:,n2),ad_H0(n2),ad_theta0(n2),ad_omega(n2),ad1_ka_drag,ad_tau_wind(:,:,n2),...
     ad_detady(:,n2),ad_dgamma(:,n2),ad_dAw(:,n2),ad_dSw(:,n2),...
     ad1_d50,ad1_d90,ad1_params] = ...
        ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Qx,ad_h(:,n2+1),bkgdn2);
    ad_ka_drag     =ad_ka_drag     +ad1_ka_drag     ;
    ad_d50         =ad_d50         +ad1_d50         ;
    ad_d90         =ad_d90         +ad1_d90         ;
    ad_params.fv   =ad_params.fv   +ad1_params.fv   ;
    ad_params.ks   =ad_params.ks   +ad1_params.ks   ;
    if(strcmp(bkgdn2.sedmodel,'vanderA'))
      ad_params.n    =ad_params.n    +ad1_params.n    ;
      ad_params.m    =ad_params.m    +ad1_params.m    ;
      ad_params.xi   =ad_params.xi   +ad1_params.xi   ;
      ad_params.alpha=ad_params.alpha+ad1_params.alpha;
    elseif(strcmp(bkgdn2.sedmodel,'dubarbier'))
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
  if(strcmp(bkgd1.sedmodel,'vanderA'))
    ad_full_params.n    (i)=ad_full_params.n    ;
    ad_full_params.m    (i)=ad_full_params.m    ;
    ad_full_params.xi   (i)=ad_full_params.xi   ;
    ad_full_params.alpha(i)=ad_full_params.alpha;
  elseif(strcmp(bkgd1.sedmodel,'dubarbier'))
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
