%
% AD symmetry test
%
clear

% % TEST: select just one i/o variable to test
% inoutvar='k';
% inoutvar='c';
% inoutvar='n ';
% inoutvar='cg';
% inoutvar='refconst';
% inoutvar='s0';
% inoutvar='gamma';
% inoutvar='theta';
% inoutvar='tharg';
% inoutvar='Hm';
% inoutvar='B';
% inoutvar='Qo';
% inoutvar='Qb';
% inoutvar='Db';
% inoutvar='Ew';
% inoutvar='Dr';
% inoutvar='Er';
% inoutvar='H';
% inoutvar='dSxydx ';
% inoutvar='Fy';
% inoutvar='Cd';
% inoutvar='urms';
% inoutvar='N';
% inoutvar='v';

load ../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define input variables
x=waves.x;
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
tauw=.001*ones(length(x),1);
detady=ones(length(x),1)*.001;
nx=length(h);

% background NL model run
[Hrms,theta,vbar,kabs,Ew,Er,Dr,bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tauw,detady);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(3*nx+4,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_h      =F(0*nx+[1:nx],i);
  tl_H0     =F(1*nx+1,i);
  tl_theta0 =F(1*nx+2,i);
  tl_omega  =F(1*nx+3,i);
  tl_ka_drag=F(1*nx+4,i);
  tl_tauw   =F(1*nx+4+[1:nx],i);
  tl_detady =F(2*nx+4+[1:nx],i);
  [tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
      tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tauw,tl_detady,bkgd);

  % AD model: g=AD*(TL*F)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauw,ad_detady] = ...
      ad_hydro_ruessink2001(tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr,bkgd);

  % create output vector
  g(0*nx+[1:nx],i)    =ad_h      ;
  g(1*nx+1,i)         =ad_H0     ;
  g(1*nx+2,i)         =ad_theta0 ;
  g(1*nx+3,i)         =ad_omega  ;
  g(1*nx+4,i)         =ad_ka_drag;
  g(1*nx+4+[1:nx],i)  =ad_tauw   ;
  g(2*nx+4+[1:nx],i)  =ad_detady ;

end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
