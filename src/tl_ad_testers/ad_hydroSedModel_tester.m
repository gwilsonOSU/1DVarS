%
% AD symmetry test
%
clear

% TEST: select just one i/o variable to test
% inoutvar='h';  %ok
% inoutvar='dh';  %ok? 10e-12
% inoutvar='hp';  %ok? 10e-12
% inoutvar='dQdx';  %ok? 10e-12
% inoutvar='qp';  %nan
% inoutvar='Qx';  %ok? 10e-12
% inoutvar='Q';  %ok? 10e-12
% inoutvar='theta';  %ok
% inoutvar='tanbeta';  %ok
% inoutvar='Hrms';  %ok
% inoutvar='kabs';  %ok
% inoutvar='Ew';  %ok
% inoutvar='Er';  %ok
% inoutvar='Dr';  %ok
% inoutvar='vbar';  %ok
% inoutvar='omega';  %ok
% inoutvar='udelta';  %ok for 1st dim and 2nd dim
% inoutvar='ubar';  %
% inoutvar='udelta_w';  %
% inoutvar='ws';  %
% inoutvar='d50';  %
% inoutvar='d90';  %
% inoutvar='k';  %
% inoutvar='detady';  %
% inoutvar='tau_wind';  %
% inoutvar='d50_8';  %
% inoutvar='ubarx';  %ok
% inoutvar='c';  %
% inoutvar='H0';  %
% inoutvar='theta0';  %
% inoutvar='ka_drag';  %
% inoutvar='dgamma';  %
% inoutvar='delta_bl';  % ok
% inoutvar='hp';


% select qtrans model to be tested
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% Use one of the duck94 test case at time step 'targetn' as background
% state.  Note, model input data are pre-cached as mat-file.
duck94Case='c'; targetn=250;
% duck94Case='b'; targetn=200;
load(['../../test_cases/duck94/obsdataCache/obsdata_case' duck94Case '.mat']);
modelinput=initModelInputs(duck94Case,grid,sedmodel);

% extract relevant model inputs for time step 'targetn'
dt=diff([hydroobs(targetn+[0 1]).dnum_est])*24*3600;
nsubsteps=1;
x=grid.x;
nx=length(x);
H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,hydroobs(targetn).dnum_est);
theta0=interp1(waves8m.dnum_est,waves8m.theta0,hydroobs(targetn).dnum_est);
omega =interp1(waves8m.dnum_est,waves8m.sigmam,hydroobs(targetn).dnum_est);
tau_wind=interp1(windEOP.dnum_est,windEOP.tau,hydroobs(targetn).dnum_est);
tau_wind=repmat(tau_wind,nx,1);
dgamma=ones(nx,1)*.01;
dAw=ones(nx,1)*01;
dSw=ones(nx,1)*01;
detady=ones(nx,1)*.0001;
tide=interp1(waves8m.dnum_est,waves8m.tide,hydroobs(targetn).dnum_est);
d50      =modelinput.d50      ;
d90      =modelinput.d90      ;
ka_drag  =modelinput.ka_drag  ;
beta0    =modelinput.beta0    ;
betaType =modelinput.betaType ;
gammaType=modelinput.gammaType;
params=modelinput.params;
h=grid.h+tide;
h(h<.75)=.75;  % min depth constraint

% OPTIONAL: override default hydro model formulations to test their TL-AD
% gammaType=2003;
% betaType='const';
% param.streamingType='v';  % choose either 'n' or 'v'

% background NL model run
[Hrms,vbar,theta,kabs,Qx,hpout,bkgd] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
                  d50,d90,params,sedmodel,gammaType,betaType,dt,nsubsteps);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(9*nx+14,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_h            =F(0*nx+[1:nx],i);
  tl_H0           =F(1*nx+1,i);
  tl_theta0       =F(1*nx+2,i);
  tl_omega        =F(1*nx+3,i);
  tl_ka_drag      =F(1*nx+4,i);
  tl_tau_wind(:,1)=F(1*nx+4+[1:nx],i);
  tl_tau_wind(:,2)=F(2*nx+4+[1:nx],i);
  tl_detady       =F(3*nx+4+[1:nx],i);
  tl_dgamma       =F(4*nx+4+[1:nx],i);
  tl_dAw          =F(5*nx+4+[1:nx],i);
  tl_dSw          =F(6*nx+4+[1:nx],i);
  tl_d50          =F(7*nx+4+[1:nx],i);
  tl_d90          =F(8*nx+4+[1:nx],i);
  tl_beta0        =F(9*nx+5,i);
  tl_params.fv    =F(9*nx+6,i);
  tl_params.ks    =F(9*nx+7,i);
  tl_params.n     =F(9*nx+8,i);
  tl_params.m     =F(9*nx+9,i);
  tl_params.xi    =F(9*nx+10,i);
  tl_params.alpha =F(9*nx+11,i);
  tl_params.lambda=F(9*nx+12,i);
  if(isfield(params,'Cc'))
    tl_params.Cc    =F(9*nx+13,i);
    tl_params.Cf    =F(9*nx+14,i);
  end

  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hpout] = ...
      tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                       tl_detady,tl_dgamma,tl_dAw,tl_dSw,...
                       tl_d50,tl_d90,tl_params,bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_beta0,ad_tau_wind,...
   ad_detady,ad_dgamma,ad_dAw,ad_dSw,...
   ad_d50,ad_d90,ad_params] = ...
      ad_hydroSedModel(tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hpout,bkgd);%,inoutvar);

  % create output vector
  g(0*nx+[1:nx],i)       =ad_h       ;
  g(1*nx+1,i)            =ad_H0      ;
  g(1*nx+2,i)            =ad_theta0  ;
  g(1*nx+3,i)            =ad_omega   ;
  g(1*nx+4,i)            =ad_ka_drag ;
  g(1*nx+4+[1:nx],i)     =ad_tau_wind(:,1);  % 1e-15
  g(2*nx+4+[1:nx],i)     =ad_tau_wind(:,2);
  g(3*nx+4+[1:nx],i)     =ad_detady  ;
  g(4*nx+4+[1:nx],i)     =ad_dgamma  ;
  g(5*nx+4+[1:nx],i)     =ad_dAw     ;
  g(6*nx+4+[1:nx],i)     =ad_dSw     ;
  g(7*nx+4+[1:nx],i)     =ad_d50     ;
  g(8*nx+4+[1:nx],i)     =ad_d90     ;
  g(9*nx+5,i)            =ad_beta0;
  g(9*nx+6,i)            =ad_params.fv;
  g(9*nx+7,i)            =ad_params.ks;
  g(9*nx+8,i)            =ad_params.n    ;
  g(9*nx+9,i)            =ad_params.m    ;
  g(9*nx+10,i)            =ad_params.xi   ;  % nan?
  g(9*nx+11,i)           =ad_params.alpha;  % nan?
  g(9*nx+12,i)           =ad_params.lambda;   % 1e-15
  if(isfield(params,'Cc'))
    g(9*nx+13,i)           =ad_params.Cc;
    g(9*nx+14,i)           =ad_params.Cf;
  end

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
