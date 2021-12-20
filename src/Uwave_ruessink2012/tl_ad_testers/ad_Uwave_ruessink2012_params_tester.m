%
% TL-AD symmetry test for Uwave_ruessink2012.m
%
clear
addpath ..
addpath ../../hydro_ruessink2001
    
load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define input variables
dt=60*60;  % one hour time step
nsubsteps=1;
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
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
sedmodel='vanderA';
params.n=1.2;
params.m=11;
params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19

% NL model run
[Hrms,theta,v,k,Ew,Er,Dr,Aw,Sw,Uw,bkgd_hydro]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma);
Hmo=1.4*Hrms;
[Aw,Sw,Uw,bkgd_uwave1]=Uwave_ruessink2012_params(Hmo,k,omega,h);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(10*nx,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n

  % % test: try with hydro model instead of Uwave
  % tl_h      =F(0*nx+0+[1:nx],i);
  % tl_H0     =F(1*nx+1       ,i);
  % tl_theta0 =F(1*nx+2       ,i);
  % tl_omega  =F(1*nx+3       ,i);
  % tl_ka_drag=F(1*nx+4       ,i);
  % tl_tauw2d(:,1) =F(1*nx+4+[1:nx],i);
  % tl_tauw2d(:,2) =F(2*nx+4+[1:nx],i);
  % tl_detady =F(3*nx+4+[1:nx],i);
  % tl_dgamma =F(4*nx+4+[1:nx],i);
  % [tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
  %     tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tauw2d,tl_detady,tl_dgamma,bkgd_hydro);
  % [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauw2d,ad_detady,ad_dgamma] = ...
  %     ad_hydro_ruessink2001(tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw,bkgd_hydro);
  % g(0*nx+0+[1:nx],i)=ad_h           ;
  % g(1*nx+1       ,i)=ad_H0          ;
  % g(1*nx+2       ,i)=ad_theta0      ;
  % g(1*nx+3       ,i)=ad_omega       ;
  % g(1*nx+4       ,i)=ad_ka_drag     ;
  % g(1*nx+4+[1:nx],i)=ad_tauw2d(:,1) ;
  % g(2*nx+4+[1:nx],i)=ad_tauw2d(:,2) ;
  % g(3*nx+4+[1:nx],i)=ad_detady      ;
  % g(4*nx+4+[1:nx],i)=ad_dgamma      ;

  tl_Hmo  =F(0*nx+[1:nx],i);
  tl_k    =F(1*nx+[1:nx],i);
  tl_omega=F(2*nx+1,i);
  tl_h    =F(2*nx+1+[1:nx],i);
  [tl_Aw,tl_Sw,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_k,tl_omega,tl_h,bkgd_uwave1);
  [ad_Hmo,ad_k,ad_omega,ad_h]=ad_Uwave_ruessink2012_params(tl_Aw,tl_Sw,tl_Uw,bkgd_uwave1);
  g(:,i)=[ad_Hmo; ad_k; ad_omega; ad_h];

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
