%
% TL-AD symmetry test for Uwave_ruessink2012.m
%
clear

% % TEST: choose an i/o variable to test
% inoutvar='kabs';
% inoutvar='A';
% inoutvar='uorb';
% inoutvar='z0';
% inoutvar='fw';
% inoutvar='Df';
% inoutvar='nubar_tb';
% inoutvar='Hm0';
% inoutvar='ht';
% inoutvar='p1';
% inoutvar='delta';
% inoutvar='phi_b';
% inoutvar='tau_wave';
% inoutvar='tau_t';
% inoutvar='nubar_tflow';
% inoutvar='abs_tau_wind';
% inoutvar='nubar_twind';
% inoutvar='nubar_twave';
% inoutvar='nubar_t';
% inoutvar='nu_tsurf';
% inoutvar='sig_s';
% inoutvar='phi_s';
% inoutvar='nums';
% inoutvar='dens';
% inoutvar='sig_b';
% inoutvar='sig_0';
% inoutvar='dsb';
% inoutvar='sgridb';
% inoutvar='dsm';
% inoutvar='sgridm';
% inoutvar='l1';
% inoutvar='l2';
% inoutvar='l3';
% inoutvar='l4';
% inoutvar='l1d';
% inoutvar='l2d';
% inoutvar='Ab';
% inoutvar='Bbp';
% inoutvar='Cbp';
% inoutvar='alphab';
% inoutvar='betab';
% inoutvar='Am';
% inoutvar='alpham';
% inoutvar='betam';
% inoutvar='alphab_bar';
% inoutvar='betab_bar';
% inoutvar='alpham_bar';
% inoutvar='betam_bar';
% inoutvar='alpha_bar';
% inoutvar='beta_bar';
% inoutvar='F';
% inoutvar='Bb';
% inoutvar='Cb';
% inoutvar='udelta';

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=20;
bkgd.ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
bkgd.ubar(2)=waves.v(i);
bkgd.k     =waves.k(i)*[cos(waves.theta(i)) sin(waves.theta(i))];
bkgd.omega =waves.sigma;
bkgd.h     =waves.h(i);
bkgd.Hrms  =waves.H(i);
bkgd.detady=.001;
bkgd.windW =[1 1]*.5;
bkgd.Dr    =waves.eps_r(i)*1030;
bkgd.fv    =.1;
bkgd.ks=.008;
bkgd.d50   =200e-6;

% bkgd NL state
[udelta,delta,udel_bkgd]=udelta_reniers2004(bkgd.ubar,...
                                            bkgd.k,...
                                            bkgd.omega,...
                                            bkgd.h,...
                                            bkgd.Hrms,...
                                            bkgd.detady,...
                                            bkgd.windW,...
                                            bkgd.Dr,...
                                            bkgd.fv,...
                                            bkgd.ks,...
                                            bkgd.d50);%,'working_nl.mat');

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(14,n);  % 1st dim is number of tl input parameters
for i=1:n

  % TL model: TL*F
  tl_ubar  =F(1:2,i)';
  tl_k     =F(3:4,i)';
  tl_omega =F(5,i);
  tl_h     =F(6,i);
  tl_Hrms  =F(7,i);
  tl_detady=F(8,i);
  tl_tau_wind=F(9:10,i)';
  tl_Dr    =F(11,i);
  tl_fv    =F(12,i);
  tl_ks    =F(13,i);
  tl_d50   =F(14,i);
  [tl_udelta,tl_delta]=tl_udelta_reniers2004(tl_ubar,tl_k,tl_omega,...
                                             tl_h,tl_Hrms,tl_detady,...
                                             tl_tau_wind,tl_Dr,tl_fv,tl_ks,tl_d50,...
                                             udel_bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_ubar,ad_k,ad_omega,ad_h,ad_Hrms,ad_detady,ad_tau_wind,ad_Dr,ad_fv,ad_ks,ad_d50] = ...
      ad_udelta_reniers2004(tl_udelta,tl_delta,udel_bkgd);%,inoutvar);
  g(:,i)=[ad_ubar(:);ad_k(:);ad_omega;ad_h;ad_Hrms;...
          ad_detady;ad_tau_wind(:);ad_Dr;ad_fv;ad_ks;ad_d50];

end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
