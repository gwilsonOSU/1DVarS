function [posterior,bkgd_final]=hydroAssimOneStep(prior,obs,dt,nsubsteps,verb)
%
% [posterior,bkgd_final]=hydroAssimOneStep(prior,obs,dt,nsubsteps,verb)
%
% Run a single time step of hydroWaveModel.m, and assimilate data to correct
% for hydrodynamic model errors.  This code (or variants of it) was formerly
% named 'assim_1dh.m' but has evolved into the current form.
%
% NOTE, this script is not usually called as a standalone function.  See
% ../hydroAssimLoop.m which uses this script in a time-loop fashion which is
% the intended use.
%
% INPUT:
%
% prior = model inputs and covariances needed for assimilation.  This needs
% to have all the fields required as input for hydroWaveModel.m (see main
% src directory for 1DVarS, which must be in your path), as well as
% background covariances needed for assimilating data.  The required
% covariances are listed below along with the variables they control:
%
%   Ch      (h, bathymetry)
%   CH0     (H0, BC rms wave height, scalar)
%   Ctheta0 (theta0, BC mean wave angle, scalar)
%   Cka     (ka_drag, bottom roughness, scalar)
%   Cgamma  (dgamma, linear correction to wave model gamma)
%   CdAw    (dAw, linear correction to wave velocity asymmetry)
%   CdSw    (dSw, linear correction to wave velocity skewness)
%
%   If you want to *not* correct one or more variables, you can set its
%   covariance to zero.  If you want *big* corrections, set a big
%   covariance.
% 
% obs = struct of observations (see prepObsData.m), with fields that can be
% any of the following obs-types:
%
%   obs.H (wave height)
%   obs.v (longshore current)
%   obs.A (wave velocity asymmetry)
%   obs.S (wave velocity skewness)
%   obs.h (total mean water depth, including tide)
%          * NOTE, assimilation of bathymetry is experimental and not recommended.
%
%   For each observation type obs.X, need sub-fields obs.X.{d,e,ind}, which
%   specify the data (d), obs-error-stdev (e), and model grid indeces (ind)
%   where the observation is located.
%
% dt,nsubsteps = time step for bathymetry update, and how many sub-steps to
% use within that time.  The substeps use identical hydrodynamic forcing,
% but reduce the time between morphodynamic updates which can improve
% stability in some cases.
%
% verb = set to 1 to show diagnostic plots while assimilating data (default 0)
%
% OUTPUT:
%
% posterior = same info as prior (input) except now updated in two ways.
% First, posterior includes 'hp' which is the bathymetry predicted at time
% t+dt.  Second, posterior includes updates to hydrodynamic forcing
% variables that were obtained by assimilating data (input 'obs').  Note the
% posterior struct omits any sub-stepping information.
%
% (OPTIONAL) bkgd_final = same as posterior, but includes all the sub-steps.
% This is sometimes necessary when running TL-AD calculations.  If you are
% using nsubsteps==1, then you don't need this.
%

if(~exist('verb'))
  verb=0;
end
if(~exist('nsubsteps'))
  nsubsteps=1;
end
nitermax=50;

% unpack some variables for convenience
Ch=prior.Ch;
Cgamma=prior.Cgamma;
CH0=prior.CH0;
Ctheta0=prior.Ctheta0;
Cka=prior.Cka;
CdAw=prior.CdAw;
CdSw=prior.CdSw;
x=prior.x;
nx=length(prior.x);

% if any obs types missing, set to empty
fld='HvhAS';
noobs.ind=[];
noobs.d=[];
noobs.e=[];
for i=1:length(fld)
  if(~isfield(obs,fld(i)))
    obs=setfield(obs,fld(i),noobs);
  end
end

% initialize bkgd hydro model state at t=0 using the prior.  Note this does
% not forecast morphology, it's just hydro at t=0
bkgd = hydroWaveModel(prior.x,prior.h,prior.H0,prior.theta0,prior.omega,prior.ka_drag,prior.beta0,prior.tau_wind,prior.detady,prior.dgamma,prior.dAw,prior.dSw,prior.gammaType,prior.betaType);
horig=prior.h;
for fld=fields(bkgd)'
  fld=cell2mat(fld);
  prior=setfield(prior,fld,getfield(bkgd,fld));
end
prior.h=horig;  % don't cut off dry points

% outer loop
eps=nan;
for n=1:nitermax

  % run t=0 TL-based prediction for this outer loop iteration, linearized
  % about most-recent background state.  The latter (bkgd) is updated witht
  % each outer loop iteration.
  fcst=prior;
  tl_h=fcst.h-bkgd.h;
  tl_H0=fcst.H0-bkgd.H0;
  tl_theta0=fcst.theta0-bkgd.theta0;
  tl_ka_drag=fcst.ka_drag-bkgd.ka_drag;
  tl_beta0  =fcst.beta0  -bkgd.beta0  ;
  tl_omega=0; % no uncertainty in omega
  tl_tau_wind=zeros(nx,2);  % no uncertainty in tauw
  tl_detady=zeros(nx,1); % no uncertainty in detady
  tl_dgamma=fcst.dgamma-bkgd.dgamma;
  tl_dAw=fcst.dAw-bkgd.dAw;
  tl_dSw=fcst.dAw-bkgd.dSw;
  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
    tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                     tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
  fcst.h=bkgd.h+tl_h;
  fcst.Hrms=bkgd.Hrms+tl_Hrms;
  fcst.vbar=bkgd.vbar+tl_vbar;
  fcst.Aw=bkgd.Aw+tl_Aw;
  fcst.Sw=bkgd.Sw+tl_Sw;

  % initialize representer sub-matrices
  %
  % LEGEND:
  %
  % R_XY = sensitivity of observation type Y, to delta-perturbations of
  % observation X.  These are sub-blocks of L*M*C*M'*L' (note, C is the
  % prior covariance).
  %
  % r_X = sensitivity of model input vector phi = [h; H0; theta0; ka_drag])
  % to delta-perturbations of observation X.  These are rows of C*M'*L'.
  %
  clear R_* r_*
  r_H=zeros(4*nx+3,length(obs.H.ind));
  r_v=zeros(4*nx+3,length(obs.v.ind));
  r_h=zeros(4*nx+3,length(obs.h.ind));
  r_A=zeros(4*nx+3,length(obs.A.ind));
  r_S=zeros(4*nx+3,length(obs.S.ind));
  R_HH=zeros(length(obs.H.ind),length(obs.H.ind));
  R_Hv=zeros(length(obs.H.ind),length(obs.v.ind));
  R_Hh=zeros(length(obs.H.ind),length(obs.h.ind));
  R_HA=zeros(length(obs.H.ind),length(obs.A.ind));
  R_HS=zeros(length(obs.H.ind),length(obs.S.ind));
  R_vH=zeros(length(obs.v.ind),length(obs.H.ind));
  R_vv=zeros(length(obs.v.ind),length(obs.v.ind));
  R_vh=zeros(length(obs.v.ind),length(obs.h.ind));
  R_vh=zeros(length(obs.v.ind),length(obs.h.ind));
  R_vA=zeros(length(obs.v.ind),length(obs.A.ind));
  R_vS=zeros(length(obs.v.ind),length(obs.S.ind));
  R_hH=zeros(length(obs.h.ind),length(obs.H.ind));
  R_hv=zeros(length(obs.h.ind),length(obs.v.ind));
  R_hh=zeros(length(obs.h.ind),length(obs.h.ind));
  R_hA=zeros(length(obs.h.ind),length(obs.A.ind));
  R_hS=zeros(length(obs.h.ind),length(obs.S.ind));
  R_AH=zeros(length(obs.A.ind),length(obs.H.ind));
  R_Av=zeros(length(obs.A.ind),length(obs.v.ind));
  R_Ah=zeros(length(obs.A.ind),length(obs.h.ind));
  R_AA=zeros(length(obs.A.ind),length(obs.A.ind));
  R_AS=zeros(length(obs.A.ind),length(obs.S.ind));
  R_SH=zeros(length(obs.S.ind),length(obs.H.ind));
  R_Sv=zeros(length(obs.S.ind),length(obs.v.ind));
  R_Sh=zeros(length(obs.S.ind),length(obs.h.ind));
  R_SA=zeros(length(obs.S.ind),length(obs.A.ind));
  R_SS=zeros(length(obs.S.ind),length(obs.S.ind));

  % compute representers for wave height: apply delta-perturbations to
  % observations of type X, to compute (a) sensitivity of model inputs
  % (r_* matrices), and (b) sensitivity of observations of type Y (R_*
  % matrices)
  for i=1:length(obs.H.ind)
    ad_Hrms =zeros(nx,1);
    ad_theta=zeros(nx,1);
    ad_vbar =zeros(nx,1);
    ad_kabs =zeros(nx,1);
    ad_Ew   =zeros(nx,1);
    ad_Er   =zeros(nx,1);
    ad_Dr   =zeros(nx,1);
    ad_Aw   =zeros(nx,1);
    ad_Sw   =zeros(nx,1);
    ad_Uw   =zeros(nx,1);
    ad_Hrms(obs.H.ind(i))=1;  % comb

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
     ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
        ad_hydroWaveModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Ew,ad_Er,ad_Dr,ad_Aw,ad_Sw,ad_Uw,bkgd);

    % multiply into covariances to get matrix of representers C*M'.  Store all
    % the variables that are uncertain.
    tl_h=Ch*ad_h;
    tl_H0=CH0*ad_H0;
    tl_theta0=Ctheta0*ad_theta0;
    tl_omega=0*ad_theta0;  % assume omega is known
    tl_ka_drag=Cka*ad_ka_drag;
    tl_tau_wind=0*ad_tau_wind;  % assume wind is known
    tl_detady=0*ad_detady;  % assume detady is known
    tl_dgamma=Cgamma*ad_dgamma;
    tl_dAw=CdAw*ad_dAw;
    tl_dSw=CdSw*ad_dSw;
    r_H(:,i)=[tl_h; tl_H0; tl_theta0; tl_ka_drag; tl_dgamma; tl_dAw; tl_dSw];

    % apply TL model and measure it to get representer matrix, M*C*M'.  Store
    % all the variables that are in the observations.
    [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
        tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                          tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
    R_HH(i,:)=tl_Hrms(obs.H.ind);
    R_Hv(i,:)=tl_vbar(obs.v.ind);
    R_Hh(i,:)=tl_h(obs.h.ind);
    R_HA(i,:)=tl_Aw(obs.A.ind);
    R_HS(i,:)=tl_Sw(obs.S.ind);

  end

  % compute representers for longshore current: same operations as in previous
  % block of code
  for i=1:length(obs.v.ind)
    ad_Hrms =zeros(nx,1);
    ad_theta=zeros(nx,1);
    ad_vbar =zeros(nx,1);
    ad_kabs =zeros(nx,1);
    ad_Ew   =zeros(nx,1);
    ad_Er   =zeros(nx,1);
    ad_Dr   =zeros(nx,1);
    ad_Aw   =zeros(nx,1);
    ad_Sw   =zeros(nx,1);
    ad_Uw   =zeros(nx,1);
    ad_vbar(obs.v.ind(i))=1;  % comb

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
     ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
        ad_hydroWaveModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Ew,ad_Er,ad_Dr,ad_Aw,ad_Sw,ad_Uw,bkgd);

    % multiply into covariances to get matrix of representers C*M'.  Store all
    % the variables that are uncertain.
    tl_h=Ch*ad_h;
    tl_H0=CH0*ad_H0;
    tl_theta0=Ctheta0*ad_theta0;
    tl_omega=0*ad_theta0;  % assume omega is known
    tl_ka_drag=Cka*ad_ka_drag;
    tl_tau_wind=0*ad_tau_wind;  % assume wind is known
    tl_detady=0*ad_detady;  % assume detady is known
    tl_dgamma=Cgamma*ad_dgamma;
    tl_dAw=CdAw*ad_dAw;
    tl_dSw=CdSw*ad_dSw;
    r_v(:,i)=[tl_h; tl_H0; tl_theta0; tl_ka_drag; tl_dgamma; tl_dAw; tl_dSw];

    % apply TL model and measure it to get representer matrix, M*C*M'.  Store
    % all the variables that are in the observations.
    [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
        tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                          tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
    R_vH(i,:)=tl_Hrms(obs.H.ind);
    R_vv(i,:)=tl_vbar(obs.v.ind);
    R_vh(i,:)=tl_h(obs.h.ind);
    R_vA(i,:)=tl_Aw(obs.A.ind);
    R_vS(i,:)=tl_Sw(obs.S.ind);

  end

  % compute representers for bathymetry: same operations as in previous block
  % of code, but bathymetry doesn't need an adjoint model run
  for i=1:length(obs.h.ind)
    ad_Hrms =zeros(nx,1);
    ad_theta=zeros(nx,1);
    ad_vbar =zeros(nx,1);
    ad_kabs =zeros(nx,1);
    ad_Ew   =zeros(nx,1);
    ad_Er   =zeros(nx,1);
    ad_Dr   =zeros(nx,1);
    ad_Aw   =zeros(nx,1);
    ad_Sw   =zeros(nx,1);
    ad_Uw   =zeros(nx,1);
    ad_h(obs.h.ind(i))=1;  % comb

    % special case: Do not need to run adjiont because all AD "inputs" are zero.
    % Because bathymetry is a model output, not a model input.

    % multiply into covariances to get matrix of representers C*M'.  Store all
    % the variables that are uncertain.
    tl_h=Ch*ad_h;
    tl_H0=CH0*ad_H0;
    tl_theta0=Ctheta0*ad_theta0;
    tl_omega=0*ad_theta0;  % assume omega is known
    tl_ka_drag=Cka*ad_ka_drag;
    tl_tau_wind=0*ad_tau_wind;  % assume wind is known
    tl_detady=0*ad_detady;  % assume detady is known
    tl_dgamma=Cgamma*ad_dgamma;
    tl_dAw=CdAw*ad_dAw;
    tl_dSw=CdSw*ad_dSw;
    r_h(:,i)=[tl_h; tl_H0; tl_theta0; tl_ka_drag; tl_dgamma; tl_dAw; tl_dSw];

    % apply TL model and measure it to get representer matrix, M*C*M'.  Store
    % all the variables that are in the observations.
    [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
        tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                          tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
    R_hH(i,:)=tl_Hrms(obs.H.ind);
    R_hv(i,:)=tl_vbar(obs.v.ind);
    R_hh(i,:)=tl_h(obs.h.ind);
    R_hA(i,:)=tl_Aw(obs.A.ind);
    R_hS(i,:)=tl_Sw(obs.S.ind);

  end

  % compute representers for wave asymmetry: same operations as in previous
  % block of code
  for i=1:length(obs.A.ind)
    ad_Hrms =zeros(nx,1);
    ad_theta=zeros(nx,1);
    ad_vbar =zeros(nx,1);
    ad_kabs =zeros(nx,1);
    ad_Ew   =zeros(nx,1);
    ad_Er   =zeros(nx,1);
    ad_Dr   =zeros(nx,1);
    ad_Aw   =zeros(nx,1);
    ad_Sw   =zeros(nx,1);
    ad_Uw   =zeros(nx,1);
    ad_Aw(obs.A.ind(i))=1;  % comb

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
     ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
        ad_hydroWaveModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Ew,ad_Er,ad_Dr,ad_Aw,ad_Sw,ad_Uw,bkgd);

    % multiply into covariances to get matrix of representers C*M'.  Store all
    % the variables that are uncertain.
    tl_h=Ch*ad_h;
    tl_H0=CH0*ad_H0;
    tl_theta0=Ctheta0*ad_theta0;
    tl_omega=0*ad_theta0;  % assume omega is known
    tl_ka_drag=Cka*ad_ka_drag;
    tl_tau_wind=0*ad_tau_wind;  % assume wind is known
    tl_detady=0*ad_detady;  % assume detady is known
    tl_dgamma=Cgamma*ad_dgamma;
    tl_dAw=CdAw*ad_dAw;
    tl_dSw=CdSw*ad_dSw;
    r_A(:,i)=[tl_h; tl_H0; tl_theta0; tl_ka_drag; tl_dgamma; tl_dAw; tl_dSw];

    % apply TL model and measure it to get representer matrix, M*C*M'.  Store
    % all the variables that are in the observations.
    [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
        tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                          tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
    R_AH(i,:)=tl_Hrms(obs.H.ind);
    R_Av(i,:)=tl_vbar(obs.v.ind);
    R_Ah(i,:)=tl_h(obs.h.ind);
    R_AA(i,:)=tl_Aw(obs.A.ind);
    R_AS(i,:)=tl_Sw(obs.S.ind);

  end

  % compute representers for wave skewness: same operations as in previous
  % block of code
  for i=1:length(obs.S.ind)
    ad_Hrms =zeros(nx,1);
    ad_theta=zeros(nx,1);
    ad_vbar =zeros(nx,1);
    ad_kabs =zeros(nx,1);
    ad_Ew   =zeros(nx,1);
    ad_Er   =zeros(nx,1);
    ad_Dr   =zeros(nx,1);
    ad_Aw   =zeros(nx,1);
    ad_Sw   =zeros(nx,1);
    ad_Uw   =zeros(nx,1);
    ad_Sw(obs.S.ind(i))=1;  % comb

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
     ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
        ad_hydroWaveModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Ew,ad_Er,ad_Dr,ad_Aw,ad_Sw,ad_Uw,bkgd);

    % multiply into covariances to get matrix of representers C*M'.  Store all
    % the variables that are uncertain.
    tl_h=Ch*ad_h;
    tl_H0=CH0*ad_H0;
    tl_theta0=Ctheta0*ad_theta0;
    tl_omega=0*ad_theta0;  % assume omega is known
    tl_ka_drag=Cka*ad_ka_drag;
    tl_tau_wind=0*ad_tau_wind;  % assume wind is known
    tl_detady=0*ad_detady;  % assume detady is known
    tl_dgamma=Cgamma*ad_dgamma;
    tl_dAw=CdAw*ad_dAw;
    tl_dSw=CdSw*ad_dSw;
    r_S(:,i)=[tl_h; tl_H0; tl_theta0; tl_ka_drag; tl_dgamma; tl_dAw; tl_dSw];

    % apply TL model and measure it to get representer matrix, M*C*M'.  Store
    % all the variables that are in the observations.
    [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
        tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                          tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);
    R_SH(i,:)=tl_Hrms(obs.H.ind);
    R_Sv(i,:)=tl_vbar(obs.v.ind);
    R_Sh(i,:)=tl_h(obs.h.ind);
    R_SA(i,:)=tl_Aw(obs.A.ind);
    R_SS(i,:)=tl_Sw(obs.S.ind);

  end

  % assemble matrices for updating
  Cd=diag([obs.H.e(:);
           obs.v.e(:);
           obs.h.e(:);
           obs.A.e(:);
           obs.S.e(:)].^2);
  CMt=[r_H r_v r_h r_A r_S];
  d=[obs.H.d(:);
     obs.v.d(:);
     obs.h.d(:);
     obs.A.d(:);
     obs.S.d(:)];
  Lu=[fcst.Hrms(obs.H.ind);
      fcst.vbar(obs.v.ind);
      fcst.h(obs.h.ind);
      fcst.Aw(obs.A.ind);
      fcst.Sw(obs.S.ind);];
  R=[R_HH R_Hv R_Hh R_HA R_HS;
     R_vH R_vv R_vh R_vA R_vS;
     R_hH R_hv R_hh R_hA R_hS;
     R_AH R_Av R_Ah R_AA R_AS;
     R_SH R_Sv R_Sh R_SA R_SS];
  obstype=[repmat('H',[length(obs.H.d) 1]);
           repmat('v',[length(obs.v.d) 1]);
           repmat('h',[length(obs.h.d) 1]);
           repmat('A',[length(obs.A.d) 1]);
           repmat('S',[length(obs.S.d) 1])];

  % compute the update.  Note, the order of parameters in the output variable
  % 'update' is determined by the order I packed them into the representers
  % (r_H, r_v, etc.)
  hedge=1; %tanh(n/5);  % reduce magnitude of update for 1st few iterations
  update=hedge*CMt*inv(R+Cd)*(d-Lu);
  % if(n>1)
  %   update=.5*(u0+update);
  % end
  % u0=update;
  posterior=prior;
  posterior.h      =prior.h      + update(0*nx+0+[1:nx]);
  posterior.H0     =prior.H0     + update(1*nx+1       );
  posterior.theta0 =prior.theta0 + update(1*nx+2       );
  posterior.ka_drag=prior.ka_drag+ update(1*nx+3       );
  posterior.dgamma =prior.dgamma + update(1*nx+3+[1:nx]);
  posterior.dAw    =prior.dAw    + update(2*nx+3+[1:nx]);
  posterior.dSw    =prior.dSw    + update(3*nx+3+[1:nx]);
  posterior.ka_drag=max(1e-4,posterior.ka_drag);  % limiter on ka_drag
  bkgd0=bkgd;

  % obtain new background state for t=0
  bkgd = hydroWaveModel(posterior.x,posterior.h,...
                        posterior.H0,posterior.theta0,posterior.omega,...
                        posterior.ka_drag,posterior.beta0,posterior.tau_wind,posterior.detady,...
                        posterior.dgamma,posterior.dAw,posterior.dSw,posterior.gammaType,posterior.betaType);

  % show the results
  if(verb)
    clf
    subplot(321), hold on
    plot(prior.x,prior.h,'r')
    plot(prior.x,prior.h-sqrt(diag(prior.Ch)),'r--')
    plot(prior.x,prior.h+sqrt(diag(prior.Ch)),'r--')
    plot(bkgd.x,bkgd.h,'b')
    set(gca,'ydir','reverse')
    ylabel('h [m]')
    % legend('prior','inverted')
    subplot(322), hold on
    plot(prior.x,prior.Hrms,'r')
    plot(bkgd.x,bkgd.Hrms,'b')
    errorbar(x(obs.H.ind),obs.H.d,obs.H.e,'ko')
    ylabel('H_{rms} [m]')
    title(['H0 prior=' num2str(prior.H0,2) 'm, final=' num2str(bkgd.H0,2) 'm']);
    subplot(323), hold on
    plot(prior.x,rad2deg(prior.theta),'r')
    plot(bkgd.x,rad2deg(bkgd.theta),'b')
    ylabel('theta [deg]')
    title(['theta0 prior=' num2str(rad2deg(prior.theta0),2) 'deg, final=' num2str(rad2deg(bkgd.theta0),2) 'deg']);
    subplot(324), hold on
    plot(prior.x,prior.vbar,'r')
    plot(bkgd.x,bkgd.vbar,'b')
    errorbar(x(obs.v.ind),obs.v.d,obs.v.e,'ko')
    ylabel('v [m/s]')
    title(['ka prior=' num2str(prior.ka_drag,3) 'm, final=' num2str(bkgd.ka_drag,3) 'm']);
    subplot(325), hold on
    plot(prior.x,prior.Aw,'r')
    plot(bkgd.x,bkgd.Aw,'b')
    errorbar(x(obs.A.ind),obs.A.d,obs.A.e,'ko')
    ylabel('Aw [-]')
    subplot(326), hold on
    plot(prior.x,prior.Sw,'r')
    plot(bkgd.x,bkgd.Sw,'b')
    errorbar(x(obs.S.ind),obs.S.d,obs.S.e,'ko')
    ylabel('Sw [-]')
    for i=1:6
      subplot(3,2,i)
      box on
      % drawAxis
    end
    pause(.1)
  end

  % check for convergence
  eps=0;
  if(trace(Ch)>0)
    eps=eps+sum((bkgd0.h-bkgd.h).^2)/trace(prior.Ch);
  end
  if(CH0>0)
    eps=eps+(bkgd0.H0-bkgd.H0)^2/prior.CH0;
  end
  if(Ctheta0>0)
    eps=eps+(bkgd0.theta0-bkgd.theta0)^2/prior.Ctheta0;
  end
  if(Cka>0)
    eps=eps+(bkgd0.ka_drag-bkgd.ka_drag)^2/prior.Cka;
  end
  if(trace(Cgamma)>0)
    eps=eps+sum((bkgd0.dgamma-bkgd.dgamma).^2)/trace(prior.Cgamma);
  end
  if(trace(CdAw)>0)
    eps=eps+sum((bkgd0.dAw-bkgd.dAw).^2)/trace(prior.CdAw);
  end
  if(trace(CdSw)>0)
    eps=eps+sum((bkgd0.dSw-bkgd.dSw).^2)/trace(prior.CdSw);
  end
  if(eps<1e-4)  % changed from 1e-3 Dec 28, 2021
    break;
  end

  disp(['iteration ' num2str(n) ', itermax = ' num2str(nitermax) ', eps = ' num2str(eps)])
end  % outer loop iterations

% finalize posterior covariances
C2=blkdiag(Ch,CH0,Ctheta0,Cka,Cgamma,CdAw,CdSw)-CMt*inv(R+Cd)*CMt';
if(min(diag(C2))<0)
  warning('C2 has negatives on diagonal!  Enforcing positive definiteness')
  C2 = forcePosDef(C2);
end
posterior.Ch      =C2(0*nx+0+[1:nx],0*nx+0+[1:nx]);
posterior.CH0     =C2(1*nx+1       ,1*nx+1       );
posterior.Ctheta0 =C2(1*nx+2       ,1*nx+2       );
posterior.Cka     =C2(1*nx+3       ,1*nx+3       );
posterior.Cgamma  =C2(1*nx+3+[1:nx],1*nx+3+[1:nx]);
posterior.CdAw    =C2(2*nx+3+[1:nx],2*nx+3+[1:nx]);
posterior.CdSw    =C2(3*nx+3+[1:nx],3*nx+3+[1:nx]);

% forecast hp for the next obs-time t+dt.  Note: the output 'posterior' only
% stores information from the final sub-step, while the 'bkgd' needs to
% include all sub-steps since it will be used in TL-AD below for propagating
% the covariance.
posterior0=posterior;
disp('running deterministic forecast')
bkgd = hydroSedModel(posterior.x,posterior.h,...
                     posterior.H0,posterior.theta0,posterior.omega,...
                     posterior.ka_drag,posterior.beta0,posterior.tau_wind,posterior.detady,...
                     posterior.dgamma,posterior.dAw,posterior.dSw,...
                     posterior.d50,posterior.d90,posterior.params,posterior.sedmodel,...
                     posterior.gammaType,posterior.betaType,...
                     dt,nsubsteps);
if(nargout==2)
  bkgd_final=bkgd;  % optional output
end
posterior = mergestruct(posterior,bkgd(end));  % update all deterministic fields, keeping only final substep
posterior.h=posterior0.h;  % exception: bkgd version has hmin cutoff
