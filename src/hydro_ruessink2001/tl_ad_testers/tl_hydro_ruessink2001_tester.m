%
% perturbation test
%
addpath ..
clear

addpath ./example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define input variables
x=waves.x;
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
tauw=ones(length(x),1)*.001;
detady=zeros(length(x),1)*.001;

% background NL model run
[Hrms,theta,vbar,kabs,Ew,Er,Dr,bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tauw,detady);

% choose perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_h      =h      *frac_tl*myrand();
tl_H0     =H0     *frac_tl*myrand();
tl_theta0 =theta0 *frac_tl*myrand();
tl_omega  =omega  *frac_tl*myrand();
tl_ka_drag=ka_drag*frac_tl*myrand();
tl_tauw   =tauw   *frac_tl*myrand();
tl_detady =detady *frac_tl*myrand();

% re-run NL model with perturbed background variables
[Hrms_prime,theta_prime,vbar_prime,kabs_prime,Ew_prime,Er_prime,Dr_prime] = ...
    hydro_ruessink2001(x,...
                       h      +tl_h      ,...
                       H0     +tl_H0     ,...
                       theta0 +tl_theta0 ,...
                       omega  +tl_omega  ,...
                       ka_drag+tl_ka_drag,...
                       tauw   +tl_tauw   ,...
                       detady +tl_detady );

% run TL model
[tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
    tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tauw,tl_detady,bkgd);

% check results
clf
vname={};
vname{end+1}='vbar';
vname{end+1}='Hrms';
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

