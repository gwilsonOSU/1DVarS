function [params,diagnostics]=bathyAssim(bkgd,bathyobs)
%
% [params,diagnostics]=bathyAssim(bkgd,bathyobs)
%
% Phase-2 assimilation.  Use bathymetry observations to correct sediment
% transport parameters.
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

% prep params_std, allow for 20% error in all parameters.  NOTE, when
% defining 'fld' below the order is important!!!
params=bkgd(1).params;
if(strcmp(bkgd(1).sedmodel,'vanderA'))
  fld={'fv','ks','n','m','xi','alpha'};
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  fld={'fv','ks','Cw','Cc','Cf','Ka'};
else
  error(['sedmodel=' num2str(bkgd(1).sedmodel) ' is not supported for full-run assimilation']);
end
for i=1:length(fld)
  params_std(i)=abs(getfield(params,fld{i}))*.2;
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
    ad_params=paramsHandler(0,bkgd(1).sedmodel,0,0,0,0,0,0);  % init ad_params struct to zero
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

    % multiply adjoint output by covariances to get the matrix of representers
    % C*M'. Note, need to do some gymnastics to convert the params struct
    % to a vector (and back). Also note, assume we only correcting for
    % unknown sediment transport params, therefore most variables are
    % simply set to zero (i.e., they have covariance=0). If correcting for
    % other things, e.g. initial bathymetry, you would have to mulitply by
    % nonzero covariance here, then measure the covariance vs. time.
    Cparam=diag(params_std.^2);
    [ad_fv,ad_ks,ad_n,ad_m,ad_xi,ad_alpha]=paramsHandler(1,bkgd(1).sedmodel,ad_params);
    if(strcmp(bkgd(1).sedmodel,'vanderA'))
      ad_pp=[ad_fv; ad_ks; ad_n; ad_m; ad_xi; ad_alpha];
    elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
      ad_pp=[ad_fv; ad_ks; ad_Cw; ad_Cc; ad_Cf; ad_Ka];
    end
    tl_pp = Cparam*ad_pp;
    tl_params=paramsHandler(0,bkgd(1).sedmodel,tl_pp(1),tl_pp(2),tl_pp(3),tl_pp(4),tl_pp(5),tl_pp(6));
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
      [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_h(:,n2+1)] = ...
          tl_hydroSedModel(tl_h(:,n2),tl_H0(n2),tl_theta0(n2),tl_omega(n2),...
                           tl_ka_drag,tl_tau_wind(:,:,n2),...
                           tl_detady(:,n2),tl_dgamma(:,n2),...
                           tl_dAw(:,n2),tl_dSw(:,n2),...
                           tl_d50,tl_d90,tl_params,bkgd(n2));
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
elseif(strcmp(bkgd(1).sedmodel,'dubarbier'))
  params.Cw = params0.Cw + update(3);
  params.Cc = params0.Cc + update(4);
  params.Cf = params0.Cf + update(5);
  params.Ka = params0.Ka + update(6);
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
