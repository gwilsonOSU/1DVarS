function [params,diagnostics]=bathyAssim(bkgd,bathyobs)
%
% [params,diagnostics]=bathyAssim(bkgd,bathyobs)
%
% Phase-2 assimilation.  Use bathymetry observations to correct sediment
% transport parameters.
%
% NOTE: This function will use the local directory /tmp/bathyAssimCache as a
% disk cache to work around memory limitation.  If the directory doesn't
% exist, it will be created, and the cached data will be deleted when
% finished..  You must have ~5GB free in /tmp, for typical runs with ~500
% time steps.
%
% INPUTS:
%
% bkgd = array of structs of model output, produced by hydroAssimLoop.m
%
% bathyobs = bathymetry observations to be assimilated, produced by prepObsData.m
%
% OUTPUTS:
%
% params = updated sediment transport parameters after assimilating data
%
% (OPTIONAL) diagnostics = matrices and such used in assimilation
%

% Phase-2 inversion code bathyAssim.m has been massaged into working with
% parfor, but is very memory-bound. Explanation: parfor allocates a copy of
% the 'bkgd' struct-array (~4.5GB) for every worker, in addition to the
% "local" copy.  I seemed to be able to get away with 2 workers on plank
% which has 62GB RAM.  This led to implementing disk caching, where the
% 'bkgd' variable is cached to disk rather than fed to parallel workers in
% RAM.  With caching, up to 12 parallel workers is fine.
parpoolN=12;
currentPool=gcp('nocreate');
if(isempty(currentPool) | currentPool.NumWorkers ~= parpoolN)
  if(~isempty(currentPool))
    delete(gcp('nocreate'));
  end
  parpool('local',parpoolN);
end
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

% prep params_std, allow for 10% error in all parameters.  NOTE, the order
% of indexes in params_std is important, it must match the ordering
% convention used by paramsHandler.m
params=bkgd(1).params;
if(strcmp(bkgd(1).sedmodel,'vanderA'))
  params_std(1)=.0101;  % default params.fv = 0.101
  params_std(2)=.00082; % default params.ks = 0.0082
  params_std(3)=.12;    % default params.n = 1.2
  params_std(4)=1.1;    % default params.m = 11
  params_std(5)=.17;    % default params.xi = 1.7
  params_std(6)=.82;    % default params.alpha = 8.2
  params_std(7)=.001;   % default params.Cc = 0.01
  params_std(8)=.003;   % default params.Cf = 0.03
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  fld={'fv','ks','Cw','Cc','Cf','Ka'};
  for i=1:length(fld)
    params_std(i)=abs(getfield(params,fld{i}))*.1;
  end
else
  error(['sedmodel=' num2str(bkgd(1).sedmodel) ' is not supported for full-run assimilation']);
end
Cd50=diag((.1*bkgd(1).d50).^2);  % optional, allow for correction to d50
Cd50=0;  % disable d50 corrections

% Init matrices needed for assimilation. Note dimensions of the full
% observation vector is no*nt x 1, where no is number of stations and nt
% is number of time points. I will break this out into a (no,nt) sized
% matrix, then later reshape it to no*nt x 1. Similar for other matrices.
% I list dimensions for each matrix below to reflect how things will later
% be reshaped in this fashion. Note I have ensured (above) that each time
% step always has the same number of observations (variable 'no')
modelnx=length(bkgd(1).x);
obsnt=length(bathyobs);  % number of observation times
modelnt=min(max([bathyobs.obsn]),length(bkgd));  % number of model time steps for a full run
obsno=length(bathyobs.measind);  % number of obs per time step
modelnp=length(params_std);
CMt =nan(modelnp,obsnt,obsno);  % np x (nt,no)
MCMt=nan(obsnt,obsno,obsnt,obsno);  % (nt,no) x (nt,no)
Mf  =nan(obsnt,obsno);  % (nt,no) x 1
d   =nan(obsnt,obsno);  % (nt,no) x 1
Cd  =nan(obsnt,obsno);  % (nt,no) x 1

% calculate representers for each observation [bathyobs.h]
for n=1:obsnt
  disp(['Computing representers for time step ' num2str(n) ' of ' num2str(obsnt)])
  tic

  clear MCMt_thisn  % parfor requires careful handling of MCMt
  parfor i=1:obsno
    if(floor(i/obsno*10)>floor((i-1)/obsno*10))
      disp(['  checkpoint ' num2str(floor(i/obsno*10)) ' of 10'])
    end

    % initialize comb for observation of h at obs-gridpoint i and obs-time n
    ad_Hrms =zeros(modelnx,1);
    ad_theta=zeros(modelnx,1);
    ad_vbar =zeros(modelnx,1);
    ad_kabs =zeros(modelnx,1);
    ad_Qx   =zeros(modelnx,1);
    ad_h   =zeros(modelnx,modelnt+1);
    ad_h(bathyobs.measind(i),bathyobs(n).obsn+1)=1;  % comb

    % initialize adjoint outputs
    % if(isempty(tmpdir))
    %   bkgd1=bkgd(1);
    % else
      bkgd1=load([tmpdir '/bkgd1.mat']);
    % end
    if(strcmp(bkgdn2.sedmodel,'vanderA'))
      ad_params=paramsHandler(0,bkgd1.sedmodel,zeros(8,1));  % init ad_params struct to zero
    elseif(strcmp(bkgdn2.sedmodel,'dubarbier'))
      ad_params=paramsHandler(0,bkgd1.sedmodel,zeros(6,1));  % init ad_params struct to zero
    end
    ad_ka_drag=0;
    ad_d50=zeros(modelnx,1);
    ad_d90=zeros(modelnx,1);
    ad_H0=zeros(modelnt,1);
    ad_theta0  =zeros(modelnt,1);
    ad_omega   =zeros(modelnt,1);
    ad_dgamma  =zeros(modelnx,modelnt);
    ad_dAw     =zeros(modelnx,modelnt);
    ad_dSw     =zeros(modelnx,modelnt);
    ad_tau_wind=zeros(modelnx,2,modelnt);
    ad_detady  =zeros(modelnx,modelnt);

    % propagate adjoint backwards from time bathyobs(n) to 1
    for n2=bathyobs(n).obsn:-1:1
      % if(isempty(tmpdir))
      %   bkgdn2=bkgd(n2);
      % else
        bkgdn2=load([tmpdir '/bkgd' num2str(n2) '.mat']);
      % end
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
        ad_params.Cc   =ad_params.Cc   +ad1_params.Cc   ;
        ad_params.Cf   =ad_params.Cf   +ad1_params.Cf   ;
      elseif(strcmp(bkgdn2.sedmodel,'dubarbier'))
        ad_params.Cw = ad_params.Cw + ad1_params.Cw;
        ad_params.Cc = ad_params.Cc + ad1_params.Cc;
        ad_params.Cf = ad_params.Cf + ad1_params.Cf;
        ad_params.Ka = ad_params.Ka + ad1_params.Ka;
      end
    end

    % multiply adjoint output by covariances to get the matrix of representers
    % C*M'. Note, need to do some gymnastics to convert the params struct
    % to a vector (and back). Also note, assume we only correcting for
    % unknown sediment transport params, therefore most variables are
    % simply set to zero (i.e., they have covariance=0). If correcting for
    % other things, e.g. initial bathymetry, you would have to mulitply by
    % nonzero covariance here, then measure the covariance vs. time.
    Cparam=diag(params_std.^2);
    if(strcmp(bkgdn2.sedmodel,'vanderA'))
      [ad_fv,ad_ks,ad_n,ad_m,ad_xi,ad_alpha,ad_Cc,ad_Cf]=paramsHandler(1,bkgdn2.sedmodel,ad_params);
      ad_pp=[ad_fv; ad_ks; ad_n; ad_m; ad_xi; ad_alpha; ad_Cc; ad_Cf];
    elseif(strcmp(bkgdn2.sedmodel,'dubarbier'))
      [ad_fv,ad_ks,ad_Cw,ad_Cc,ad_Cf,ad_Ka]=paramsHandler(1,bkgdn2.sedmodel,ad_params);
      ad_pp=[ad_fv; ad_ks; ad_Cw; ad_Cc; ad_Cf; ad_Ka];
    end
    tl_pp = Cparam*ad_pp;
    tl_params=paramsHandler(0,bkgdn2.sedmodel,tl_pp);
    tl_h = zeros(modelnx,modelnt+1);
    tl_h(:,1)      =0*ad_h(:,1);  % bathymetry at time n=1
    tl_H0          =0*ad_H0      ;
    tl_theta0      =0*ad_theta0  ;
    tl_omega       =0*ad_omega   ;
    tl_ka_drag     =0*ad_ka_drag ;
    tl_tau_wind    =0*ad_tau_wind;
    tl_detady      =0*ad_detady  ;
    tl_dgamma      =0*ad_dgamma  ;
    tl_dAw         =0*ad_dAw     ;
    tl_dSw         =0*ad_dSw     ;
    tl_d50         =Cd50*ad_d50  ;
    tl_d90         =0*ad_d90     ;

    % since tl_params is constant, we need only measure it once for the
    % observation we're considering
    CMt(:,n,i)=tl_pp;

    % Run TL model forwards in time for full simulation period
    for n2=1:modelnt
      % if(isempty(tmpdir))
      %   bkgdn2=bkgd(n2);
      % else
        bkgdn2=load([tmpdir '/bkgd' num2str(n2) '.mat']);
      % end
      [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_h(:,n2+1)] = ...
          tl_hydroSedModel(tl_h(:,n2),tl_H0(n2),tl_theta0(n2),tl_omega(n2),...
                           tl_ka_drag,tl_tau_wind(:,:,n2),...
                           tl_detady(:,n2),tl_dgamma(:,n2),...
                           tl_dAw(:,n2),tl_dSw(:,n2),...
                           tl_d50,tl_d90,tl_params,bkgdn2);
    end  % TL forward time loop (index n2)

    % measure the TL output to get representers for observation (n,i)
    MCMt_thisn{i} = tl_h(bathyobs.measind,[bathyobs.obsn]);

  end  % loop (index i) over observations at time n

  % weird variable MCMt_thisn{i}(n,i) was needed to get parfor to run.  Unpack
  % it now.
  for i=1:obsno
    MCMt(n,i,:,:)=MCMt_thisn{i}';
  end

  % store observations for obs time n
  Mf(n,:)=bkgd(bathyobs(n).obsn).hp(bathyobs(n).h.ind);  % (nt,no) x 1
  d (n,:)=bathyobs(n).h.d;  % (nt,no) x 1
  Cd(n,:)=bathyobs(n).h.e.^2;  % (nt,no) x 1

  disp(['  step ' num2str(n) ' finished in ' num2str(toc) 'sec'])
end  % bathyobs time loop (index n)

% reshape matrices s.t. observations are treated as one big vector
CMt =reshape(CMt ,[modelnp obsnt*obsno]);  % np x (nt,no)
MCMt=reshape(MCMt,[obsnt*obsno obsnt*obsno]);  % (nt,no) x (nt,no)
Mf  =reshape(Mf  ,[obsnt*obsno 1]);  % (nt,no) x 1
d   =reshape(d   ,[obsnt*obsno 1]);  % (nt,no) x 1
Cd  =reshape(Cd  ,[obsnt*obsno 1]);  % (nt,no) x 1

% assimilation step, update the parameters
ind=find(~isnan(d));
update=CMt(:,ind)*inv(MCMt(ind,ind)+diag(Cd(ind)))*(d(ind)-Mf(ind));
params0=bkgd(1).params;
params.fv   =params0.fv   +update(1);
params.ks   =params0.ks   +update(2);
if(strcmp(bkgd(1).sedmodel,'vanderA'))
  params.n    =params0.n    +update(3);
  params.m    =params0.m    +update(4);
  params.xi   =params0.xi   +update(5);
  params.alpha=params0.alpha+update(6);
  params.Cc   =params0.Cc   +update(7);
  params.Cf   =params0.Cf   +update(8);
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  params.Cw = params0.Cw + update(3);
  params.Cc = params0.Cc + update(4);
  params.Cf = params0.Cf + update(5);
  params.Ka = params0.Ka + update(6);
end

% clean up disk cache
disp('Deleting cached bkgd struct data')
nt=length(bkgd);
for n=1:nt
  unix(['rm ' tmpdir '/bkgd' num2str(n) '.mat']);
end

% pack output struct with assimilation diagnostics
if(nargout==2)
  diagnostics.CMt    =CMt    ;
  diagnostics.MCMt   =MCMt   ;
  diagnostics.Mf     =Mf     ;
  diagnostics.d      =d      ;
  diagnostics.Cd     =Cd     ;
  diagnostics.update =update ;
  diagnostics.params0=params0;
  diagnostics.params_std=params_std;
end
