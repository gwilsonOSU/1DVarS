%
% perturbation test of tl_qtrans_vanderA.m
%
clear

% % TEST-CODE: choose a test variable.  These are listed in the order they
% % appear in the TL code
% inoutvar='Cd';
% inoutvar='uE2';
% inoutvar='ucr';
% inoutvar='urms';
% inoutvar='Hmo';
% inoutvar='aw';
% inoutvar='Ur';
% inoutvar='psi';
% inoutvar='earg';
% inoutvar='dens';
% inoutvar='Sk';
% inoutvar='As';
% inoutvar='VW';
% inoutvar='uAV';
% inoutvar='Dstar';
% inoutvar='Asb';
% inoutvar='Ass';
% inoutvar='Ufact';
% inoutvar='slopeFact';
% inoutvar='C';
% inoutvar='Dh';
% inoutvar='qtot';

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
x=waves.x(i);
h=waves.h(i);
Hrms=waves.H(i);
ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
ubar(2)=waves.v(i);
kabs=waves.k(i);
omega=waves.sigma;
detady=.001;
windW =[1 1]*.5;
Dr=waves.eps_r(i)*1030;
fv=1;
d50=180e-6;
d90=400e-6;
ws=ws_brownLawler(d50);
param.alphab=1.6;
param.facua =1.0;
tanbeta=.01;  % scalar
theta=waves.theta(i);

[q,bkgd] = qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,kabs,...
                                 omega,theta,ubar,Dr,param);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F=eps*rand(13,n);  % 1st dim is number of tl input parameters
for i=1:n

  % TL model: TL*F
  tl_d50         =F(1,i);
  tl_d90         =F(2,i);
  tl_h           =F(3,i);
  tl_tanbeta     =F(4,i);
  tl_Hrms        =F(5,i);
  tl_kabs        =F(6,i);
  tl_omega       =F(7,i);
  tl_theta       =F(8,i);
  tl_ubar        =F(9:10,i)';
  tl_Dr          =F(11,i);
  tl_param.facua =F(12,i);
  tl_param.alphab=F(13,i);
  tl_qs = tl_qtrans_soulsbyVanRijn(tl_d50,...
                                   tl_d90,...
                                   tl_h,...
                                   tl_tanbeta,...
                                   tl_Hrms,...
                                   tl_kabs,...
                                   tl_omega,...
                                   tl_theta,...
                                   tl_ubar,...
                                   tl_Dr,...
                                   tl_param,...
                                   bkgd); %,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_d50,ad_d90,ad_h,ad_tanbeta,...
   ad_Hrms,ad_kabs,ad_omega,ad_theta,ad_ubar,...
   ad_Dr,ad_param] = ad_qtrans_soulsbyVanRijn(tl_qs,bkgd); %,inoutvar);
  g(:,i)=[ad_d50;ad_d90;ad_h;ad_tanbeta;...
          ad_Hrms;ad_kabs;ad_omega;ad_theta;ad_ubar(:);...
          ad_Dr;ad_param.facua;ad_param.alphab];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
