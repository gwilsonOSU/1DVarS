%
% perturbation test
%
clear

% TEST: specify output variable to override
% outvarlist={};
% outvarlist{end+1}='Hmo';
% outvarlist{end+1}='Aw';
% outvarlist{end+1}='Sw';
% outvarlist{end+1}='Uw';
% outvarlist{end+1}='k(:,1)';
% outvarlist{end+1}='k(:,2)';
% outvarlist{end+1}='c ';
% outvarlist{end+1}='ubarx ';
% outvarlist{end+1}='vbar';
% outvarlist{end+1}='ws ';
% outvarlist{end+1}='xb ';
% outvarlist{end+1}='ubar(:,1) ';
% outvarlist{end+1}='ubar(:,2) ';
% outvarlist{end+1}='delta_bl';
% outvarlist{end+1}='udelta(:,1)';
% outvarlist{end+1}='udelta(:,2)';
% outvarlist{end+1}='udelta_w(:,1)';
% outvarlist{end+1}='udelta_w(:,2)';
% outvarlist{end+1}='tanbeta ';
% outvarlist{end+1}='ur(:,2)';
% outvarlist{end+1}='udelta_w(:,1)';
% outvarlist{end+1}='udelta_w(:,2)';
% outvarlist{end+1}='Q0 ';
% outvarlist{end+1}='d50';
% outvarlist{end+1}='d90';
% outvarlist{end+1}='h';
% outvarlist{end+1}='tanbeta';
% outvarlist{end+1}='Hrms';
% outvarlist{end+1}='kabs';
% outvarlist{end+1}='omega';
% outvarlist{end+1}='delta_bl';
% outvarlist{end+1}='Q1';
% outvarlist{end+1}='Q ';
% outvarlist{end+1}='Qx ';
% outvarlist{end+1}='dh ';
% outvarlist{end+1}='dQdx ';
% outvarlist{end+1}='dh ';
% outvarlist{end+1}='qp ';
% outvarlist{end+1}='dh ';
% outvarlist{end+1}='hp ';
% for ii=1:length(outvarlist)
%   outvar=outvarlist{ii};
%   disp(outvar)

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
                  d50,d90,params,sedmodel,gammaType,betaType,dt,nsubsteps);%,outvar);

% choose perturbations
frac_tl = 1e-4;
myrand=@(sz)2*(rand(sz)-.5);
tl_h            =h            *frac_tl.*myrand(size(h            ));
tl_H0           =H0           *frac_tl.*myrand(size(H0           ));
tl_theta0       =theta0       *frac_tl.*myrand(size(theta0       ));
tl_omega        =omega        *frac_tl.*myrand(size(omega        ));
tl_ka_drag      =ka_drag      *frac_tl.*myrand(size(ka_drag      ));
tl_beta0        =beta0        *frac_tl.*myrand(size(beta0        ));
tl_tau_wind     =tau_wind     *frac_tl.*myrand(size(tau_wind     ));
tl_detady       =detady       *frac_tl.*myrand(size(detady       ));
tl_dgamma       =dgamma       *frac_tl.*myrand(size(dgamma       ));
tl_dAw          =dAw          *frac_tl.*myrand(size(dAw          ));
tl_dSw          =dSw          *frac_tl.*myrand(size(dSw          ));
tl_d50          =d50          *frac_tl.*myrand(size(d50          ));
tl_d90          =d90          *frac_tl.*myrand(size(d90          ));
tl_params.fv    =params.fv    *frac_tl*myrand(1);
tl_params.ks    =params.ks    *frac_tl*myrand(1);
tl_params.lambda=params.lambda*frac_tl*myrand(1);
if(strcmp(sedmodel,'dubarbier'))
  tl_params.Cw  =params.Cw    *frac_tl*myrand(1);
  tl_params.Cc  =params.Cc    *frac_tl*myrand(1);
  tl_params.Cf  =params.Cf    *frac_tl*myrand(1);
  tl_params.Ka  =params.Ka    *frac_tl*myrand(1);
elseif(strcmp(sedmodel,'vanderA'))
  tl_params.n    =params.n    *frac_tl*myrand(1);
  tl_params.m    =params.m    *frac_tl*myrand(1);
  tl_params.xi   =params.xi   *frac_tl*myrand(1);
  tl_params.alpha=params.alpha*frac_tl*myrand(1);
  if(isfield(params,'Cc'))
    tl_params.Cc =params.Cc   *frac_tl*myrand(1);
    tl_params.Cf =params.Cf   *frac_tl*myrand(1);
  end
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  tl_params.alphab=params.alphab*frac_tl*myrand(1);
  tl_params.facua=params.facua*frac_tl*myrand(1);
end

% filter the spatially variable noise
m=20;
tl_h            =filtfilt(ones(m,1)/m,1,tl_h            );
tl_tau_wind(:,1)=filtfilt(ones(m,1)/m,1,tl_tau_wind(:,1));
tl_tau_wind(:,2)=filtfilt(ones(m,1)/m,1,tl_tau_wind(:,2));
tl_detady       =filtfilt(ones(m,1)/m,1,tl_detady       );
tl_dgamma       =filtfilt(ones(m,1)/m,1,tl_dgamma       );
tl_dAw          =filtfilt(ones(m,1)/m,1,tl_dAw          );
tl_dSw          =filtfilt(ones(m,1)/m,1,tl_dSw          );
tl_d50          =filtfilt(ones(m,1)/m,1,tl_d50          );
tl_d90          =filtfilt(ones(m,1)/m,1,tl_d90          );

% pack perturbed variables into params struct for running NL model
params1.fv=params.fv+tl_params.fv;
params1.ks=params.ks+tl_params.ks;
params1.lambda=params.lambda+tl_params.lambda;
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
  if(isfield(params,'Cc'))
    params1.Cc   =params.Cc   +tl_params.Cc   ;
    params1.Cf   =params.Cf   +tl_params.Cf   ;
  end
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

% drawnow;
% pause(.1);
% end
