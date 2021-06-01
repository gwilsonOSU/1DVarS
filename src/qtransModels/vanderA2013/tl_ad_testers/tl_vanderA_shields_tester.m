% NOTE: at time of writing, had not yet completed qtrans_vanderA.m.  So I
% just copied all the code up to the point of calling vanderA_shields.m, and
% took it from there.
%
clear

physicalConstants;

load ../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
ubar(2)=waves.v(i);
k     =waves.k(i)*[cos(waves.theta(i)) sin(waves.theta(i))];
omega =waves.sigma;
h     =waves.h(i);
Hrms  =waves.H(i);
detady=.001;
windW =[1 1]*.5;
Dr    =waves.eps_r(i)*1030;
fv    =.1;
d50   =200e-6;

% other user inputs required for qtrans_vanderA.m
d50=200e-6;
d90=300e-6;
kabs=sqrt(k(1)^2+k(2)^2);
ws=ws_brownLawler(d50);
udelta=udelta_reniers2004(ubar,k,omega,h,Hrms,detady,windW,Dr,fv,d50);

%------------------------------------
% begin code copied from qtrans_vanderA.m
%------------------------------------

physicalConstants;

% fixed constants
delta=0.2;  % vanderA's approximation, see Fig. 14 and discussion thereof
nt=1000;  % for calculation of intra-wave velocity

% derived params
Hmo=Hrms*1.4;

% intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% relevant to van Der A
t=linspace(0,omega/2/pi,nt);
kabs=sqrt(k(1)^2+k(2)^2);
[uw,wksp]=Uwave_ruessink2012(omega*t,Hmo,kabs,omega,h);
r_r2012  =wksp.r  ;
phi_r2012=wksp.phi;
uhat=sqrt(.5*mean(uw.^2));   % rms wave velocity for full wave cycle, eqn 8
uhatc=max(uw);  % max velocity under crest, see text below eqn 7
uhatt=-min(uw);  % max velocity under trough (positive valued)

% timing of wave velocity direction change, crest, and trough, based on
% Ruessink et al 2012.
Tc=asin(-r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)))/omega;   % duration of crest
Tt=2*pi/omega-Tc;     % duration of trough
T=Tc+Tt;  % full wave period
[~,icu_guess]=max(uw);  % 1st guess
[~,itu_guess]=min(uw);  % 1st guess
Tcu=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(icu_guess));
TTt=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(itu_guess));
Ttu=TTt-Tc;  % duration of deceleration under trough

% critical shields param
Dstar=(g*(s-1)/nu^2)^(1/3)*d50;
theta_cr=.3/(1+1.2*Dstar)+0.055*(1-exp(-.02*Dstar));

% constant parameter mu (eqn A2)
if(d50<=.15e-3)
  mu=6;
elseif(.15e-3<d50&d50<.2e-3)
  mu=6-5*(d50-.15e-3)/(.05e-3);
else
  mu=1;
end

% wave velocity moments
ahat=uhat.*T/(2*pi);  % eqn 9
utildecr=.5*sqrt(2)*uhatc;  % eqn 10
utildetr=.5*sqrt(2)*uhatt;  % eqn 11

% ripples, Appendix B
if(d50<.22e-3)
  meta=0.55;
  mlambda=.73;
elseif(.22e-3<=d50&d50<0.3e-3)
  meta=.55+.45*(d50-.22e-3)/(.3e-3-.22e-3);
  mlambda=.73+.27*(d50-.22e-3)/(.3e-3-.22e-3);
else
  meta=1;
  mlambda=1;
end
psihatc=(1.27*uhatc)^2/((s-1)*g*d50);  % 1.27 is for irregular flow
psihatt=(1.27*uhatt)^2/((s-1)*g*d50);
psihat=max(psihatc,psihatt);
if(psihat<=190)
  neta=1;
elseif(190<psihat&psihat<240)
  neta=.5*(1+cos(pi*(psihat-190)/(240-190)));
else
  neta=0;
end
nlambda=neta;
eta=ahat*meta*neta*(.275-.022*psihat^.42);
lambda=ahat*mlambda*nlambda*(1.97-0.44*psihat^.21);

%-------------------------------
% end of copied-code from qtrans_vanderA.m.  Now proceed with TL testing
%-------------------------------

% base NL solution
udeltaabs=sqrt(udelta(1)^2+udelta(2)^2);
[theta_av,ksd,ksw,fd,fw,...
          branch_A1,branch_A4,branch_A5] = ...
    vanderA_shields(d50,d90,udeltaabs,uhat,delta,mu,eta,lambda,ahat);

% TL-perturbations
tl_fac=.001;
for n=1:100

  % define perturbations
  myrand=@().5*(1-rand(1));
  tl_d50   =d50      *myrand()*tl_fac;
  tl_d90   =d90      *myrand()*tl_fac;
  tl_udelta=udeltaabs*myrand()*tl_fac;
  tl_delta =delta    *myrand()*tl_fac;
  tl_uhat  =uhat     *myrand()*tl_fac;
  tl_mu    =mu       *myrand()*tl_fac;
  tl_eta   =eta      *myrand()*tl_fac;
  tl_lambda=lambda   *myrand()*tl_fac;
  tl_ahat  =ahat     *myrand()*tl_fac;

  % NL perturbed
  [theta_av1,ksd1,ksw1,fd1,fw1,...
          branch_A11,branch_A41,branch_A51] = ...
      vanderA_shields(d50   +tl_d50   ,...
                      d90   +tl_d90   ,...
                      udeltaabs+tl_udelta,...
                      uhat  +tl_uhat  ,...
                      delta,...
                      mu    +tl_mu    ,...
                      eta   +tl_eta   ,...
                      lambda+tl_lambda,...
                      ahat  +tl_ahat  );
  tl_theta_av_true(n)=theta_av1-theta_av;
  tl_ksd_true(n)     =ksd1     -ksd     ;
  tl_ksw_true(n)     =ksw1     -ksw     ;
  tl_fd_true(n)      =fd1      -fd      ;
  tl_fw_true(n)      =fw1      -fw      ;

  % TL
  [tl_theta_av(n),tl_ksd(n),tl_ksw(n),tl_fd(n),tl_fw(n),A] = ...
      tl_vanderA_shields(tl_d50,tl_d90,tl_udelta,tl_uhat,tl_mu,tl_eta,tl_lambda,tl_ahat,...
                         d50,d90,udeltaabs,uhat,delta,mu,eta,lambda,ahat,...
                         ksd,ksw,fd,fw,...
                         branch_A1,branch_A4,branch_A5);

end


clf
subplot(3,2,1)
plot(tl_theta_av_true,tl_theta_av,'.')
title('theta av');
subplot(3,2,2)
plot(tl_ksd_true     ,tl_ksd     ,'.')
title('ksd');
subplot(3,2,3)
plot(tl_ksw_true     ,tl_ksw     ,'.')
title('ksw');
subplot(3,2,4)
plot(tl_fd_true      ,tl_fd      ,'.')
title('fd');
subplot(3,2,5)
plot(tl_fw_true      ,tl_fw      ,'.')
title('fw');
for i=1:5
  subplot(3,2,i)
  axis equal tight
  ax=max(abs(axis));
  axis([-1 1 -1 1]*ax)
  hold on
  plot([0 ax],[0 ax],'k--')
end
