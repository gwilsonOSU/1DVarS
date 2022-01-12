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

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% select dubarbier or vanderA qtrans model
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

if(strcmp(sedmodel,'dubarbier'))
  params.Cw=0.00483 ;
  params.Cc=0.02002 ;
  params.Cf=0.01173 ;
  params.Ka=0.631e-4;
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;
  params.m=11;
  params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
  params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
  % params.Cc=0.01;  % comment out to use option of automatic above-WBL transport (recommended)
  % params.Cf=0.01;
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% define input variables
dt=60*60;  % one hour time step
nsubsteps=1;
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
params.ks=0.0082;  % undertow bed roughness calibration parameter
params.lambda=1.5;
x=waves.x;
nx=length(x);
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
h(h<=.5)=.5;
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
tau_wind=ones(nx,2)*.0001;
detady=ones(nx,1)*.0001;
dgamma=zeros(nx,1);
dAw   =ones(nx,1)*+.01;
dSw   =ones(nx,1)*-.01;
d50=180e-6*ones(nx,1);
d90=400e-6*ones(nx,1);
gammaType=2003;
betaType='const';
beta0=0.1;

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
