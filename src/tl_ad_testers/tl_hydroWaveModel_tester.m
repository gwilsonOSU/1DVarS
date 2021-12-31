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

% define input variables
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
beta0=0.1;
gammaType=2003;
betaType='const';

% background NL model run
[Hrms,vbar,theta,kabs,Ew,Er,Dr,Aw,Sw,Uw,bkgd] = ...
    hydroWaveModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,gammaType,betaType);

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

% re-run NL model with perturbed background variables
[Hrms_prime,vbar_prime,theta_prime,kabs_prime,Ew_prime,Er_prime,Dr_prime,Aw_prime,Sw_prime,Uw_prime] = ...
    hydroWaveModel(x,...
                   h        +tl_h        ,...
                   H0       +tl_H0       ,...
                   theta0   +tl_theta0   ,...
                   omega    +tl_omega    ,...
                   ka_drag  +tl_ka_drag  ,...
                   beta0    +tl_beta0    ,...
                   tau_wind +tl_tau_wind ,...
                   detady   +tl_detady   ,...
                   dgamma   +tl_dgamma   ,...
                   dAw      +tl_dAw      ,...
                   dSw      +tl_dSw      ,...
                   gammaType,...
                   betaType);

% run TL model
[tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
    tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_beta0,tl_tau_wind,...
                     tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);%,outvar)

% check results
clf
vname={};
vname{end+1}='Hrms';
vname{end+1}='vbar';
vname{end+1}='theta';
vname{end+1}='kabs';
vname{end+1}='Ew';
vname{end+1}='Er';
vname{end+1}='Dr';
vname{end+1}='Aw';
vname{end+1}='Sw';
vname{end+1}='Uw';
for i=1:length(vname)
  subplot(length(vname)/2,2,i)
  hold on
  plot(eval(['tl_' vname{i}]))
  plot(eval([vname{i} '_prime - ' vname{i}]))
  title(vname{i})
  if(i==1)
    legend('TL predicted','NL diff')
  end
end

