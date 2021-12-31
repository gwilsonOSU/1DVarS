%
% perturbation test
%
clear

% TEST: specify output variable to override
% outvar='k';
% outvar='omega';
% outvar='h';
% outvar='Hrms';
% outvar='detady';
% outvar='tau_wind(:,1)';
% outvar='Dr';  %ok
% outvar='params.fv';
% outvar='params.ks';
% outvar='d50';
% outvar='vbar';  %ok
% outvar='ubarx';  %ok
% outvar='hp';  %ok!
% outvar='h';  %ok
% outvar='tanbeta';  %ok
% outvar='Hrms';  %ok
% outvar='kabs';  %ok
% outvar='omega';  %input, ok
% outvar='udelta_w(:,1)';  % biased, but ok is due to udelta_reniers
% outvar='udelta_w(:,2)';  %ok
% outvar='ws';
% outvar='Aw';  %ok
% outvar='Sw';  %ok
% outvar='Uw';  %ok

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
  params.Cc=0.01;
  params.Cf=0.03;
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% define input variables
dt=60*60;  % one hour time step
nsubsteps=1;
params.fv=0.1;
params.ks=0.0082;
x=waves.x;
nx=length(x);
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
tau_wind=ones(nx,2)*0;
detady=zeros(nx,1);
dgamma=ones(nx,1)*+.01;
dAw   =ones(nx,1)*+.01;
dSw   =ones(nx,1)*-.01;
d50=180e-6*ones(nx,1);
d90=400e-6*ones(nx,1);
beta0=0.1;
gammaType=2003;
betaType='const';

% background NL model run
[Hrms,vbar,theta,kabs,Qx,hpout,bkgd] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
                  d50,d90,params,sedmodel,gammaType,betaType,dt,nsubsteps);%,outvar);

% choose perturbations
frac_tl = 0.001;
myrand=@()2*(rand(1)-.5);
tl_h      =h      *frac_tl*myrand();
tl_H0     =H0     *frac_tl*myrand();
tl_theta0 =theta0 *frac_tl*myrand();
tl_omega  =omega  *frac_tl*myrand();
tl_ka_drag=ka_drag*frac_tl*myrand();
tl_beta0  =beta0  *frac_tl*myrand();
tl_tau_wind  =  tau_wind*frac_tl*myrand();
tl_detady =detady *frac_tl*myrand();
tl_dgamma =dgamma *frac_tl*myrand();
tl_dAw    =dAw    *frac_tl*myrand();
tl_dSw    =dSw    *frac_tl*myrand();
tl_d50    =d50    *frac_tl*myrand();
tl_d90    =d90    *frac_tl*myrand();
tl_params.fv  =params.fv   *frac_tl*myrand();
tl_params.ks  =params.ks   *frac_tl*myrand();
if(strcmp(sedmodel,'dubarbier'))
  tl_params.Cw  =params.Cw*frac_tl*myrand();
  tl_params.Cc  =params.Cc*frac_tl*myrand();
  tl_params.Cf  =params.Cf*frac_tl*myrand();
  tl_params.Ka  =params.Ka*frac_tl*myrand();
elseif(strcmp(sedmodel,'vanderA'))
  tl_params.n    =params.n    *frac_tl*myrand();
  tl_params.m    =params.m    *frac_tl*myrand();
  tl_params.xi   =params.xi   *frac_tl*myrand();
  tl_params.alpha=params.alpha*frac_tl*myrand();
  tl_params.Cc   =params.Cc   *frac_tl*myrand();
  tl_params.Cf   =params.Cf   *frac_tl*myrand();
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  tl_params.alphab=params.alphab*frac_tl*myrand();
  tl_params.facua=params.facua*frac_tl*myrand();
end

% pack perturbed variables into params struct for running NL model
params1.fv=params.fv+tl_params.fv;
params1.ks=params.ks+tl_params.ks;
if(strcmp(sedmodel,'dubarbier'))
  params1.Cw  =params.Cw+tl_params.Cw;
  params1.Cc  =params.Cc+tl_params.Cc;
  params1.Cf  =params.Cf+tl_params.Cf;
  params1.Ka  =params.Ka+tl_params.Ka;
elseif(strcmp(sedmodel,'vanderA'))
  params1.n    =params.n    +tl_params.n    ;
  params1.m    =params.m    +tl_params.m    ;
  params1.xi   =params.xi   +tl_params.xi   ;
  params1.alpha=params.alpha+tl_params.alpha;
  params1.Cc   =params.Cc   +tl_params.Cc   ;
  params1.Cf   =params.Cf   +tl_params.Cf   ;
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params1.alphab=params.alphab+tl_params.alphab;
  params1.facua =params.facua +tl_params.facua ;
end

% re-run NL model with perturbed background variables
[Hrms_prime,vbar_prime,theta_prime,kabs_prime,Qx_prime,hpout_prime] = ...
    hydroSedModel(x,...
                  h       +tl_h       ,...
                  H0      +tl_H0      ,...
                  theta0  +tl_theta0  ,...
                  omega   +tl_omega   ,...
                  ka_drag +tl_ka_drag ,...
                  beta0   +tl_beta0   ,...
                  tau_wind+tl_tau_wind,...
                  detady  +tl_detady  ,...
                  dgamma  +tl_dgamma  ,...
                  dAw     +tl_dAw     ,...
                  dSw     +tl_dSw     ,...
                  d50     +tl_d50     ,...
                  d90     +tl_d90     ,...
                  params1,...
                  sedmodel,gammaType,betaType,dt,nsubsteps);%,outvar);

% run TL model
[tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Qx,tl_hpout] = ...
    tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                     tl_detady,tl_dgamma,tl_dAw,tl_dSw,...
                     tl_d50,tl_d90,tl_params,bkgd);%,outvar);

% check results
clf
vname={};
vname{end+1}='Hrms';
vname{end+1}='vbar';
vname{end+1}='theta';
vname{end+1}='kabs';
vname{end+1}='Qx';
vname{end+1}='hpout';
for i=1:length(vname)
  subplot(length(vname),1,i)
  hold on
  plot(eval(['tl_' vname{i}]))
  plot(eval([vname{i} '_prime - ' vname{i}]))
  title(vname{i})
  if(i==1)
    legend('TL predicted','NL diff')
  end
end

