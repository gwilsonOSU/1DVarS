%
% AD symmetry test
%
clear

% TEST: select just one i/o variable to test
% inoutvar='h';

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define input variables
dt=60*60;  % one hour time step
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

% background NL model run
[Hrms,vbar,theta,kabs,Ew,Er,Dr,Aw,Sw,Uw,bkgd] = ...
    hydroWaveModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,dAw,dSw);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(9*nx+9,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_h       =F(0*nx+[1:nx],i);
  tl_H0      =F(1*nx+1,i);
  tl_theta0  =F(1*nx+2,i);
  tl_omega   =F(1*nx+3,i);
  tl_ka_drag =F(1*nx+4,i);
  tl_tau_wind(:,1)=F(1*nx+4+[1:nx],i);
  tl_tau_wind(:,2)=F(2*nx+4+[1:nx],i);
  tl_detady  =F(3*nx+4+[1:nx],i);
  tl_dgamma  =F(4*nx+4+[1:nx],i);
  tl_dAw     =F(5*nx+4+[1:nx],i);
  tl_dSw     =F(6*nx+4+[1:nx],i);

  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
      tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,...
                       tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd);

  % AD model: g=AD*(TL*F)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tau_wind,...
   ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
      ad_hydroWaveModel(tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw,bkgd);

  % create output vector
  g(0*nx+[1:nx],i)       =ad_h       ;
  g(1*nx+1,i)            =ad_H0      ;
  g(1*nx+2,i)            =ad_theta0  ;
  g(1*nx+3,i)            =ad_omega   ;
  g(1*nx+4,i)            =ad_ka_drag ;
  g(1*nx+4+[1:nx],i)=ad_tau_wind(:,1);
  g(2*nx+4+[1:nx],i)=ad_tau_wind(:,2);
  g(3*nx+4+[1:nx],i)     =ad_detady  ;
  g(4*nx+4+[1:nx],i)     =ad_dgamma  ;
  g(5*nx+4+[1:nx],i)     =ad_dAw     ;
  g(6*nx+4+[1:nx],i)     =ad_dSw     ;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
