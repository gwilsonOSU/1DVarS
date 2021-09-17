%
% perturbation test of tl_udelta_reniers2004.m
%
clear

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
bkgd.ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
bkgd.ubar(2)=waves.v(i);
bkgd.k     =waves.k(i)*[cos(waves.theta(i)) sin(waves.theta(i))];
bkgd.omega =waves.sigma;
bkgd.h     =waves.h(i);
bkgd.Hrms  =waves.H(i);
bkgd.detady=.001;
bkgd.tau_wind =[1 1]*.001;
bkgd.Dr    =waves.eps_r(i)*1030;
bkgd.fv    =.1;
bkgd.d50   =200e-6;

% bkgd NL state
[udelta,~,udel_bkgd]=udelta_reniers2004(bkgd.ubar,...
                                      bkgd.k,...
                                      bkgd.omega,...
                                      bkgd.h,...
                                      bkgd.Hrms,...
                                      bkgd.detady,...
                                      bkgd.tau_wind,...
                                      bkgd.Dr,...
                                      bkgd.fv,...
                                      bkgd.d50);%,'working_nl.mat');

verb=0;
for niter=1:100

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
tl_tau_wind=(bkgd.tau_wind -bkgd.tau_wind *(frac_tl*myrand()+1));
tl_Dr    =(bkgd.Dr    -bkgd.Dr    *(frac_tl*myrand()+1));
tl_fv    =(bkgd.fv    -bkgd.fv    *(frac_tl*myrand()+1));

% compute perturbation
udelta_prime=udelta_reniers2004(bkgd.ubar  +tl_ubar  ,...
                                bkgd.k     +tl_k     ,...
                                bkgd.omega +tl_omega ,...
                                bkgd.h     +tl_h     ,...
                                bkgd.Hrms  +tl_Hrms  ,...
                                bkgd.detady+tl_detady,...
                                bkgd.tau_wind +tl_tau_wind ,...
                                bkgd.Dr    +tl_Dr    ,...
                                bkgd.fv    +tl_fv    ,...
                                bkgd.d50   +tl_d50   );
tl_udelta(:,niter)=tl_udelta_reniers2004(tl_ubar,tl_k,tl_omega,tl_h,tl_Hrms,tl_detady,...
                                tl_tau_wind,tl_Dr,tl_fv,tl_d50,udel_bkgd);
tl_udelta_true(:,niter)=udelta_prime-udelta;

end

% comparison plot
clf
plot(tl_udelta_true,tl_udelta,'.')
hold on
axis equal tight
ax=axis;
plot(ax([1 2]),ax([1 2]),'k--')
hold off
title(['udelta comparison (should be 1-1)'])
xlabel('true [m/s]')
ylabel('predicted [m/s]')

return;

% test code: compare each tl variable
nl0=load('working_nl.mat');
nl1=load('working_nlp.mat');
tl=load('working_tl.mat');
vname={};
vname{end+1}='kabs';
vname{end+1}='A';
vname{end+1}='uorb';
vname{end+1}='ks';
vname{end+1}='z0';
vname{end+1}='fw';
vname{end+1}='Df';
vname{end+1}='nubar_tb';
vname{end+1}='Hm0';
vname{end+1}='ht';
vname{end+1}='p1';
vname{end+1}='delta';
vname{end+1}='phi_b';
vname{end+1}='tau_wind';
vname{end+1}='tau_wave';
vname{end+1}='tau_t';
vname{end+1}='nubar_tflow ';
vname{end+1}='abs_tau_wind ';
vname{end+1}='nubar_twind ';
vname{end+1}='nubar_twave ';
vname{end+1}='nubar_t ';
vname{end+1}='nu_tsurf ';
vname{end+1}='sig_s ';
vname{end+1}='phi_s ';
vname{end+1}='sig_b ';
vname{end+1}='sig_0';
vname{end+1}='Ab';
vname{end+1}='Bbp';
vname{end+1}='Cbp';
vname{end+1}='Am';
vname{end+1}='alphab_bar';
vname{end+1}='betab_bar';
vname{end+1}='alpham_bar';
vname{end+1}='betam_bar';
vname{end+1}='alpha_bar';
vname{end+1}='beta_bar';
vname{end+1}='F';
vname{end+1}='Bb';
vname{end+1}='Cb';
vname{end+1}='udelta';
for i=1:length(vname)
  disp(vname{i})
  disp(['tl: ' num2str(getfield(tl,['tl_' vname{i}]))]);
  disp(['nl: ' num2str(getfield(nl1,vname{i}) ...
                       - getfield(nl0,vname{i}))]);
  pause;
end
