%
% perturbation test of tl_udelta_reniers2004.m
%
clear

% TEST: choose an i/o variable to test
outvar='kabs';  %ok
outvar='A';  %ok
outvar='uorb';  %ok
outvar='ks';  %ok
outvar='z0';  %ok
outvar='fw';  %ok
% outvar='nubar_tb';  %ok
% outvar='Hm0';  %ok
% outvar='ht';  %ok
% outvar='p1';  %ok
% outvar='delta';  %ok
% outvar='phi_b';  %ok
% outvar='w1';
% outvar='tau_wave';
% outvar='nubar_twave';
% outvar='sig_s';
% outvar='phi_s';
% outvar='nums';
% outvar='dens';
% outvar='sig_b';  %ok
% outvar='sig_0';
% outvar='l1d';  %ok
% outvar='l2d';  %ok
% outvar='l1';
% outvar='l2';
% outvar='l3';
% outvar='l4';
% outvar='Ab';  %ok
% outvar='Bbp';
% outvar='Cbp';
% outvar='alphab';
% outvar='betab';
% outvar='Am';
% outvar='alpham';
% outvar='betam';
% outvar='alphab_bar';
% outvar='betab_bar';
% outvar='alpham_bar';
% outvar='betam_bar';
% outvar='alpha_bar';
% outvar='beta_bar';
% outvar='tau_t';  %ok
% outvar='F';  %ok
% outvar='Df';  %ok
% outvar='k';  %ok
% outvar='omega';  %ok
% outvar='Bb';  %some scatter but obviously correct
% outvar='Cb';  %ok
% outvar='udelta';

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
                                            bkgd.d50);%,outvar);

for iter=1:100

  % choose reasonable perturbations
  frac_tl = 0.01;
  myrand=@()2*(rand(1)-.5);
  tl_ubar  =(bkgd.ubar  -bkgd.ubar  *(frac_tl*myrand()+1));
  tl_k     =(bkgd.k     -bkgd.k     *(frac_tl*myrand()+1));
  tl_omega =(bkgd.omega -bkgd.omega *(frac_tl*myrand()+1));
  tl_d50   =(bkgd.d50   -bkgd.d50   *(frac_tl*myrand()+1));
  tl_h     =(bkgd.h     -bkgd.h     *(frac_tl*myrand()+1));
  tl_Hrms  =(bkgd.Hrms  -bkgd.Hrms  *(frac_tl*myrand()+1));
  tl_detady=(bkgd.detady-bkgd.detady*(frac_tl*myrand()+1));
  tl_windW =(bkgd.windW -bkgd.windW *(frac_tl*myrand()+1));
  tl_Dr    =(bkgd.Dr    -bkgd.Dr    *(frac_tl*myrand()+1));
  tl_fv    =(bkgd.fv    -bkgd.fv    *(frac_tl*myrand()+1));
  tl_ks    =(bkgd.ks    -bkgd.ks    *(frac_tl*myrand()+1));

  % compute perturbation
  udelta_prime=udelta_reniers2004(bkgd.ubar  +tl_ubar  ,...
                                  bkgd.k     +tl_k     ,...
                                  bkgd.omega +tl_omega ,...
                                  bkgd.h     +tl_h     ,...
                                  bkgd.Hrms  +tl_Hrms  ,...
                                  bkgd.detady+tl_detady,...
                                  bkgd.windW +tl_windW ,...
                                  bkgd.Dr    +tl_Dr    ,...
                                  bkgd.fv    +tl_fv    ,...
                                  bkgd.ks    +tl_ks    ,...
                                  bkgd.d50   +tl_d50   );%,outvar);
  tl_udelta_true(:,iter)=udelta_prime-udelta;
  tl_udelta(:,iter)=tl_udelta_reniers2004(tl_ubar,tl_k,tl_omega,tl_h,tl_Hrms,tl_detady,tl_windW,tl_Dr,tl_fv,tl_ks,tl_d50,udel_bkgd);%,outvar);

end

% comparison plot, scatter
clf
plot(tl_udelta_true',tl_udelta','.')
hold on
axis equal tight
ax=axis;
plot(ax([1 2]),ax([1 2]),'k--')
hold off
title(['udelta comparison (should be 1-1)'])
xlabel('true [m/s]')
ylabel('predicted [m/s]')
