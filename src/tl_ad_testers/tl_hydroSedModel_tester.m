%
% perturbation test
%
clear

load ../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
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
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% define input variables
params.dt=60*60;  % one hour time step
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
x=waves.x;
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
windW=ones(length(x),2)*2;
detady=zeros(length(x),1)*.001;
d50=180e-6;
d90=400e-6;

% background NL model run
[Hrms,vbar,theta,kabs,Q,hp,bkgd]=hydroSedModel(x,h,H0,theta0,omega,ka_drag,...
                                               windW,detady,d50,d90,params,sedmodel);

% choose perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_h      =h      *frac_tl*myrand();
tl_H0     =H0     *frac_tl*myrand();
tl_theta0 =theta0 *frac_tl*myrand();
tl_omega  =omega  *frac_tl*myrand();
tl_ka_drag=ka_drag*frac_tl*myrand();
tl_windW  =windW  *frac_tl*myrand();
tl_detady =detady *frac_tl*myrand();
tl_d50    =d50    *frac_tl*myrand();
tl_d90    =d90    *frac_tl*myrand();
tl_params.fv  =params.fv   *frac_tl*myrand();
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
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  tl_params.alphab=params.alphab*frac_tl*myrand();
  tl_params.facua=params.facua*frac_tl*myrand();
end

% re-run NL model with perturbed background variables
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
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params1.alphab=params.alphab+tl_params.alphab;
  params1.facua =params.facua +tl_params.facua ;
end
params1.fv=params.fv+tl_params.fv;
params1.dt=params.dt;
[Hrms_prime,vbar_prime,theta_prime,kabs_prime,Q_prime] = ...
    hydroSedModel(x,...
                  h      +tl_h      ,...
                  H0     +tl_H0     ,...
                  theta0 +tl_theta0 ,...
                  omega  +tl_omega  ,...
                  ka_drag+tl_ka_drag,...
                  windW  +tl_windW  ,...
                  detady +tl_detady ,...
                  d50    +tl_d50    ,...
                  d90    +tl_d90    ,...
                  params1,...
                  sedmodel);

% run TL model
[tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp] = ...
    tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,...
                     tl_windW,tl_detady,tl_d50,tl_d90,tl_params,bkgd);

% check results
clf
vname={};
vname{end+1}='vbar';
vname{end+1}='Q';
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

