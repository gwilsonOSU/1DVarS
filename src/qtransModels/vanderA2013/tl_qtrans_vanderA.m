function [tl_qs,tl_all]=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_delta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd)%,outvar)
%
% TL code for qtrans_vanderA.m
%

% this wrapper loop serves to handle vector inputs
nx=length(tl_h);
for i=1:nx
  if(bkgd(i).Hmo==0)  % ignore masked points
    tl_qs(i)=0;
  else
    if(nargout==2)
      [tl_qs(i),tl_all(i)]= ...
          tl_qtrans_vanderA_main(tl_d50(i),tl_d90(i),tl_h(i),tl_tanbeta(i),tl_Hrms(i),tl_kabs(i),...
                                 tl_omega,tl_udelta(i,:),tl_delta(i),tl_ws(i),tl_Aw(i),tl_Sw(i),tl_Uw(i),tl_param,bkgd(i));%,outvar);
    else
      tl_qs(i)= ...
          tl_qtrans_vanderA_main(tl_d50(i),tl_d90(i),tl_h(i),tl_tanbeta(i),tl_Hrms(i),tl_kabs(i),...
                                 tl_omega,tl_udelta(i,:),tl_delta(i),tl_ws(i),tl_Aw(i),tl_Sw(i),tl_Uw(i),tl_param,bkgd(i));%,outvar);
    end
  end
end
tl_qs=tl_qs(:);

end  % end of wrapper function, start of main function

function [tl_qs,tl_all]=tl_qtrans_vanderA_main(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_delta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd)%,outvar)

physicalConstants;

% % v1, break out NL background vars. This version uses eval() and is a bit slow.
% fld=fields(bkgd);
% for i=1:length(fld)
%   eval([fld{i} ' = bkgd.' fld{i} ';']);
% end
% alpha=bkgd.alpha;
% mu=bkgd.mu;

% v2, Break out NL background vars. Note this hard-coded version runs
% significantly faster than if I use eval() to dynamically load all bkgd
% variables
Dstar       =bkgd.Dstar       ;
Hmo         =bkgd.Hmo         ;
Hrms        =bkgd.Hrms        ;
Omegac      =bkgd.Omegac      ;
Omegacc     =bkgd.Omegacc     ;
Omegact     =bkgd.Omegact     ;
Omegat      =bkgd.Omegat      ;
Omegatc     =bkgd.Omegatc     ;
Omegatt     =bkgd.Omegatt     ;
Pc          =bkgd.Pc          ;
Pt          =bkgd.Pt          ;
T           =bkgd.T           ;
Tc          =bkgd.Tc          ;
Tcu         =bkgd.Tcu         ;
Tt          =bkgd.Tt          ;
Ttu         =bkgd.Ttu         ;
ahat        =bkgd.ahat        ;
alpha       =bkgd.alpha       ;
alphaw      =bkgd.alphaw      ;
branch_A1   =bkgd.branch_A1   ;
branch_A4   =bkgd.branch_A4   ;
branch_A5   =bkgd.branch_A5   ;
c           =bkgd.c           ;
d50         =bkgd.d50         ;
d90         =bkgd.d90         ;
deltasc     =bkgd.deltasc     ;
deltast     =bkgd.deltast     ;
eta         =bkgd.eta         ;
etawc       =bkgd.etawc       ;
etawt       =bkgd.etawt       ;
fd          =bkgd.fd          ;
fw          =bkgd.fw          ;
fwc         =bkgd.fwc         ;
fwdc        =bkgd.fwdc        ;
fwdt        =bkgd.fwdt        ;
fwt         =bkgd.fwt         ;
tauwRe      =bkgd.tauwRe      ;
h           =bkgd.h           ;
kabs        =bkgd.kabs        ;
ksd         =bkgd.ksd         ;
ksw         =bkgd.ksw         ;
lambda      =bkgd.lambda      ;
meta        =bkgd.meta        ;
mlambda     =bkgd.mlambda     ;
mu          =bkgd.mu          ;
neta        =bkgd.neta        ;
nlambda     =bkgd.nlambda     ;
nt          =bkgd.nt          ;
omega       =bkgd.omega       ;
param       =bkgd.param       ;
phi_r2012   =bkgd.phi_r2012   ;
r_r2012     =bkgd.r_r2012     ;
psed        =bkgd.psed        ;
psihat      =bkgd.psihat      ;
psihatc     =bkgd.psihatc     ;
psihatt     =bkgd.psihatt     ;
qs          =bkgd.qs          ;
qsc         =bkgd.qsc         ;
qst         =bkgd.qst         ;
t           =bkgd.t           ;
theta_av    =bkgd.theta_av    ;
theta_cr    =bkgd.theta_cr    ;
thetac      =bkgd.thetac      ;
thetacx     =bkgd.thetacx     ;
thetatx     =bkgd.thetatx     ;
thetahatc   =bkgd.thetahatc   ;
thetahatt   =bkgd.thetahatt   ;
thetat      =bkgd.thetat      ;
ucrabs      =bkgd.ucrabs      ;
ucrvec      =bkgd.ucrvec      ;
udabs       =bkgd.udabs       ;
udelta      =bkgd.udelta      ;
delta       =bkgd.delta       ;
uhat        =bkgd.uhat        ;
uhatc       =bkgd.uhatc       ;
uhatt       =bkgd.uhatt       ;
utildecr    =bkgd.utildecr    ;
utildetr    =bkgd.utildetr    ;
utrabs      =bkgd.utrabs      ;
utrvec      =bkgd.utrvec      ;
uw          =bkgd.uw          ;
ws          =bkgd.ws          ;
wsc         =bkgd.wsc         ;
wst         =bkgd.wst         ;
uwave_wksp  =bkgd.uwave_wksp  ;
c1          =bkgd.c1          ;
Aw          =bkgd.Aw          ;
Sw          =bkgd.Sw          ;
b           =bkgd.b           ;
RR          =bkgd.RR          ;
worb1c      =bkgd.worb1c      ;
worb1t      =bkgd.worb1t      ;
worb2c      =bkgd.worb2c      ;
worb2t      =bkgd.worb2t      ;
t1ca        =bkgd.t1ca        ;
t1ta        =bkgd.t1ta        ;
t1cb        =bkgd.t1cb        ;
t1tb        =bkgd.t1tb        ;
t1c         =bkgd.t1c         ;
t1t         =bkgd.t1t         ;
t2cb        =bkgd.t2cb        ;
t2tb        =bkgd.t2tb        ;
t2ca        =bkgd.t2ca        ;
t2ta        =bkgd.t2ta        ;
t2c         =bkgd.t2c         ;
t2t         =bkgd.t2t         ;
worbc       =bkgd.worbc       ;
worbt       =bkgd.worbt       ;
phiuc       =bkgd.phiuc       ;
phidc       =bkgd.phidc       ;
icu_guess   =bkgd.icu_guess   ;
itu_guess   =bkgd.itu_guess   ;
if(~isfield(param,'streamingType') || param.streamingType=='v')
  fwd=bkgd.fwd;
elseif(param.streamingType=='n')
  f25    =bkgd.f25    ;
  theta25=bkgd.theta25;
  r      =bkgd.r      ;
  fws    =bkgd.fws    ;
end
if(~isfield(param,'nosusp') || param.nosusp==0)
  eps_s        =bkgd.eps_s          ;
  eps_b        =bkgd.eps_b          ;
  uwmo         =bkgd.uwmo           ;
  tanbeta      =bkgd.tanbeta        ;
end
if(isfield(param,'Cc'))  % option-1 for above-WBL transport
  qs2          =bkgd.qs2            ;
  qs3          =bkgd.qs3            ;
elseif(~isfield(param,'nosusp') || param.nosusp==0)  % option-2 for above-WBL transport
  tan_phi      =bkgd.tan_phi        ;
  slopeFact    =bkgd.slopeFact      ;
  Lt     =bkgd.Lt      ;
  Lc     =bkgd.Lc      ;
  wfracc =bkgd.wfracc  ;
  wfract =bkgd.wfract  ;
  utot   =bkgd.utot    ;
  Omegai =bkgd.Omegai  ;
  Omega  =bkgd.Omega   ;
  Ksdenom1=bkgd.Ksdenom1 ;
  Kbdenom1=bkgd.Kbdenom1 ;
  Ksdenom2=bkgd.Ksdenom2 ;
  Kbdenom2=bkgd.Kbdenom2 ;
  Ki      =bkgd.Ki       ;
  K       =bkgd.K        ;
end
qsVdA        =bkgd.qsVdA          ;
qsCc         =bkgd.qsCc           ;
qsCf         =bkgd.qsCf           ;

% derived params
tl_Hmo=tl_Hrms*1.4;
tl_c=-omega/kabs^2*tl_kabs + tl_omega/kabs;

% note, van der A specifies to use significant orbital velocity amplitude,
% while Ruessink et al. (2012) uses rms.  Uwave_ruessink2012() follows the
% Ruessink et al. (2012) convention, so I need to revert to van der A's
% convention here as a special case.
tl_Uw=1.4*tl_Uw;

% intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% relevant to van Der A
[tl_uw,tl_r_r2012,tl_phi_r2012]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,omega*t,uwave_wksp);
uw2mean=mean(uw.^2);
tl_uw2mean=mean(2*uw.*tl_uw);
% uhat=sqrt(2*mean(uw.^2));   % rms wave velocity for full wave cycle, eqn 8
tl_uhat=.5./sqrt(2*uw2mean)*2.*tl_uw2mean;

% timing of wave velocity direction change, crest, and trough, based on
% Ruessink et al 2012.
asinarg = -r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2));
tl_asinarg = -tl_r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)) ...
    - r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))*tl_phi_r2012 ...
    + r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2))^2*( .5/sqrt(1-r_r2012^2)*(-2*r_r2012*tl_r_r2012) );
% phiuc=asin(asinarg)/omega;   % phase at first upcrossing
tl_phiuc = 1./sqrt(1-asinarg^2)/omega*tl_asinarg - asin(asinarg)/omega^2*tl_omega;
% phidc=pi-phiuc;  % phase at first downcrossing
tl_phidc=-tl_phiuc;
% tuc=phiuc/omega;  % time of first upcrossing
tl_tuc = tl_phiuc/omega - phiuc/omega^2*tl_omega;
% tdc=phidc/omega;  % time of first downcrossing
tl_tdc = tl_phidc/omega - phidc/omega^2*tl_omega;
% tcr=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(icu_guess));  % time of crest
tl_tcr = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,omega,r_r2012,phi_r2012,t(icu_guess));
% ttr=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(itu_guess));  % time of trough
tl_ttr = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,omega,r_r2012,phi_r2012,t(itu_guess));
% T=2*pi/omega;  % full wave period
tl_T=-2*pi/omega^2*tl_omega;
% Tc = tdc-tuc;  % duration of crest
tl_Tc = tl_tdc-tl_tuc;
% Tt=T-Tc;     % duration of trough
tl_Tt=tl_T-tl_Tc;
% Ttu=ttr-tdc;  % duration of deceleration under trough
tl_Ttu=tl_ttr-tl_tdc;
% Tcu=tcr-tuc;  % duration of acceleration under crest
tl_Tcu=tl_tcr-tl_tuc;

% for crest velocities, best I can do without being too fancy and screwing
% up the AD code is just to consider the TL velocity at the location of the
% NL model crest/trough
[~,ic]=max(uw);
[~,it]=min(uw);
tl_uhatc = tl_uw(ic);
tl_uhatt = -tl_uw(it);

% critical shields param, Soulsby
tl_Dstar=(g*(s-1)/nu^2)^(1/3)*tl_d50;
% theta_cr=.3/(1+1.2*Dstar) + 0.055*( 1 - exp(-.02*Dstar) );
tl_theta_cr = -.3/(1+1.2*Dstar)^2*(1.2*tl_Dstar) ...
    - 0.055*exp(-.02*Dstar)*(-.02*tl_Dstar);

% constant parameter mu (eqn A2)
if(d50<=.15e-3)
  tl_mu=0;
elseif(.15e-3<d50&d50<.2e-3)
  tl_mu=-5*tl_d50/.05e-3;
else
  tl_mu=0;
end

% wave velocity moments
tl_ahat = tl_uhat.*T/(2*pi) ...
          + uhat.*tl_T/(2*pi);
tl_utildecr=.5*sqrt(2)*tl_uhatc;
tl_utildetr=.5*sqrt(2)*tl_uhatt;

% ripples, Appendix B
if(d50<.22e-3)
  tl_meta=0;
  tl_mlambda=0;
elseif(.22e-3<=d50&d50<0.3e-3)
  tl_meta=.45*tl_d50/(.3e-3-.22e-3);
  tl_mlambda=.27*tl_d50/(.3e-3-.22e-3);
else
  tl_meta=0;
  tl_mlambda=0;
end
tl_psihatc = 1.27^2*2*uhatc*tl_uhatc/((s-1)*g*d50) ...
    - (1.27*uhatc)^2/((s-1)*g*d50^2)*tl_d50;
tl_psihatt = 1.27^2*2*uhatt*tl_uhatt/((s-1)*g*d50) ...
    - (1.27*uhatt)^2/((s-1)*g*d50^2)*tl_d50;
if(psihatc>psihatt)
  tl_psihat=tl_psihatc;
else
  tl_psihat=tl_psihatt;
end
if(psihat<=190)
  tl_neta=0;
elseif(190<psihat&psihat<240)
  tl_neta = -.5*sin(pi*(psihat-190)/(240-190)) ...
            *(pi*tl_psihat/(240-190));
else
  tl_neta=0;
end
tl_nlambda=tl_neta;
tl_eta = tl_ahat*meta*neta*(.275-.022*psihat^.42) ...
         + ahat*tl_meta*neta*(.275-.022*psihat^.42) ...
         + ahat*meta*tl_neta*(.275-.022*psihat^.42) ...
         - ahat*meta*neta*.022*.42*psihat^(.42-1)*tl_psihat;
tl_lambda = tl_ahat*mlambda*nlambda*(1.97-0.44*psihat^.21) ...
    + ahat*tl_mlambda*nlambda*(1.97-0.44*psihat^.21) ...
    + ahat*mlambda*tl_nlambda*(1.97-0.44*psihat^.21) ...
    - ahat*mlambda*nlambda*0.44*.21*psihat^(.21-1)*tl_psihat;

% shields parameter related parameters.  Requires solving a 5-eqn nonlinear
% system, so this is done in its own code
if(sqrt(udelta(1)^2+udelta(2)^2)==0)
  tl_udabs=0;
else
  tl_udabs=.5/sqrt(udelta(1)^2+udelta(2)^2)*( ...
      2*udelta(1)*tl_udelta(1) ...
      + 2*udelta(2)*tl_udelta(2) );
end
[tl_theta_av,tl_ksd,tl_ksw,tl_fd,tl_fw] = ...
    tl_vanderA_shields(tl_d50,tl_d90,tl_udabs,tl_uhat,tl_delta,...
                       tl_mu,tl_eta,tl_lambda,tl_ahat,...
                       d50,d90,udabs,uhat,delta,mu,eta,lambda,ahat,...
                       ksd,ksw,fd,fw,...
                       branch_A1,branch_A4,branch_A5);

% other BBL derived parameters
tl_alpha = tl_udabs./(udabs+uhat) ...
    - udabs./(udabs+uhat)^2*( tl_udabs + tl_uhat );
c1=2.6;  % cited in text as van der A et al. (2011)
if(ahat/ksw>1.587)  % eqn 21

  argc2=2*Tcu/Tc;
  tl_argc2 = 2*tl_Tcu/Tc - 2*Tcu/Tc^2*tl_Tc;
  argc1=argc2^c1*ahat/ksw;
  tl_argc1 = c1*argc2^(c1-1)*ahat/ksw*tl_argc2 ...
      + argc2^c1*tl_ahat/ksw ...
      - argc2^c1*ahat/ksw^2*tl_ksw;
  argc=5.21*argc1^(-.19);
  tl_argc = -.19*5.21*argc1^(-1.19)*tl_argc1;
  % fwc=.00251*exp(argc);
  tl_fwc=.00251*exp(argc)*tl_argc;

  argt2=2*Ttu/Tt;
  tl_argt2 = 2*tl_Ttu/Tt - 2*Ttu/Tt^2*tl_Tt;
  argt1=argt2^c1*ahat/ksw;
  tl_argt1 = c1*argt2^(c1-1)*ahat/ksw*tl_argt2 ...
      + argt2^c1*tl_ahat/ksw ...
      - argt2^c1*ahat/ksw^2*tl_ksw;
  argt=5.21*argt1^(-.19);
  tl_argt = -.19*5.21*argt1^(-1.19)*tl_argt1;
  % fwt=.00251*exp(argt);
  tl_fwt=.00251*exp(argt)*tl_argt;

else
  tl_fwc=0;
  tl_fwt=0;
end
tl_fwdc = tl_alpha*fd ...
          + alpha*tl_fd ...
          - tl_alpha*fwc ...
          + (1-alpha)*tl_fwc;
tl_fwdt = tl_alpha*fd ...
       + alpha*tl_fd ...
       - tl_alpha*fwt ...
       + (1-alpha)*tl_fwt;
tl_ucrvec = tl_utildecr*[+1 0] + tl_udelta;
tl_utrvec = tl_utildetr*[-1 0] + tl_udelta;
tl_ucrabs = .5/sqrt(ucrvec(1)^2+ucrvec(2)^2)*( ...
    2*ucrvec(1)*tl_ucrvec(1) ...
    + 2*ucrvec(2)*tl_ucrvec(2) );
tl_utrabs = .5/sqrt(utrvec(1)^2+utrvec(2)^2)*( ...
    2*utrvec(1)*tl_utrvec(1) ...
    + 2*utrvec(2)*tl_utrvec(2) );
tl_thetac = .5/((s-1)*g)*( tl_fwdc.*ucrabs.^2/d50 ...
                           + 2*fwdc.*ucrabs*tl_ucrabs/d50 ...
                           - fwdc.*ucrabs.^2/d50^2*tl_d50 );
tl_thetat = .5/((s-1)*g)*( tl_fwdt.*utrabs.^2/d50 ...
                           + 2*fwdt.*utrabs*tl_utrabs/d50 ...
                           - fwdt.*utrabs.^2/d50^2*tl_d50 );  % ad symmetric

% boundary layer streaming, either with VDA13 model or Nielsen
alphaw = 4/(3*pi);
if(~isfield(param,'streamingType') || param.streamingType=='v')
  % fwd = alpha*fd+(1-alpha)*fw;
  tl_fwd = tl_alpha*fd ...
         + alpha*tl_fd ...
         - tl_alpha*fw ...
         + (1-alpha)*tl_fw;
  % tauwRe = fwd*alphaw*uhat^3/2/c;
  tl_tauwRe = .5*alphaw*( ...
      tl_fwd*uhat^3/c ...
      + 3*fwd*uhat^2/c*tl_uhat ...
      - fwd*uhat^3/c^2*tl_c );
elseif(param.streamingType=='n')
  % f25 = exp(5.5*(2.5*d50/ahat)^.2-6.3);
  tl_f25 = exp(5.5*(2.5*d50/ahat)^2-6.3)*2*5.5*(2.5*d50/ahat)*( ...
      2.5*tl_d50/ahat - 2.5*d50/ahat^2*tl_ahat );
  % theta25 = 0.5*f25*(ahat*omega)^2/((s-1)*g*d50);
  tl_theta25 = 0.5*tl_f25*(ahat*omega)^2/((s-1)*g*d50) ...
      + 2*0.5*f25*ahat*(omega)^2/((s-1)*g*d50)*tl_ahat ...
      + 2*0.5*f25*(ahat)^2*omega/((s-1)*g*d50)*tl_omega ...
      - 0.5*f25*(ahat*omega)^2/((s-1)*g*d50^2)*tl_d50;
  % r = 170*sqrt(max(0,theta25-0.05))*d50;
  if(theta25-0.05<=0)
    tl_r=0;
  else
    tl_r = 170*sqrt(theta25-0.05)*tl_d50 ...
         + .5*170/sqrt(theta25-0.05)*tl_theta25;
  end
  if(lambda>0)
    tl_r = tl_r + 2*8*eta/lambda*tl_eta - 8*eta^2/lambda^2*tl_lambda;
  end
  % fws = exp(5.5*(r/ahat)^.2-6.3);
  tl_fws = exp(5.5*(r/ahat)^.2-6.3)*2*5.5*(r/ahat)*( tl_r/ahat - r/ahat.^2.*tl_ahat );
  % tauwRe = fws*alphaw*uhat^3/2/c;
  tl_tauwRe = tl_fws*alphaw*uhat^3/2/c ...
      + 3*fws*alphaw*uhat^2/2/c*tl_uhat ...
      - fws*alphaw*uhat^3/2/c^2*tl_c;
elseif(param.streamingType=='0')
  tl_tauwRe=0;
else
  error('must provide param.streamingType as either ''n'' or ''v''')
end

% apply streaming to bottom stress
% streamingEffect = tauwRe/((s-1)*g*d50);  % eqns 15 and 22
tl_streamingEffect = 1/((s-1)*g)*( tl_tauwRe/d50 - tauwRe/d50^2*tl_d50 );
tl_thetacx = tl_thetac*ucrvec(1)/ucrabs ...
    + thetac*tl_ucrvec(1)/ucrabs ...
    - thetac*ucrvec(1)/ucrabs^2*tl_ucrabs ...
    + tl_streamingEffect;
tl_thetatx = tl_thetat*utrvec(1)/utrabs ...
    + thetat*tl_utrvec(1)/utrabs ...
    - thetat*utrvec(1)/utrabs^2*tl_utrabs ...
    + tl_streamingEffect;

% sheet flow layer thickness, Appendix C
tl_thetahatc = .5/((s-1)*g)*( tl_fwdc.*uhatc.^2/d50 ...
                           + 2*fwdc.*uhatc*tl_uhatc/d50 ...
                           - fwdc.*uhatc.^2/d50^2*tl_d50 );
tl_thetahatt = .5/((s-1)*g)*( tl_fwdt.*uhatt.^2/d50 ...
                           + 2*fwdt.*uhatt*tl_uhatt/d50 ...
                           - fwdt.*uhatt.^2/d50^2*tl_d50 );
if(d50<=.15e-3)
  tl_deltasc=tl_d50*25*thetahatc ...
      + d50*25*tl_thetahatc;
  tl_deltast=tl_d50*25*thetahatt ...
      + d50*25*tl_thetahatt;
elseif(.15e-3<d50&d50<.2e-3)
  tl_deltasc = tl_d50*(25-12*(d50-.15e-3)/.05e-3) ...
      - d50*12*tl_d50/.05e-3;
  tl_deltast = tl_d50*(25-12*(d50-.15e-3)/.05e-3) ...
      - d50*12*tl_d50/.05e-3;
else
  tl_deltasc = tl_d50*13*thetahatc ...
      + d50*13*tl_thetahatc;
  tl_deltast = tl_d50*13*thetahatt ...
      + d50*13*tl_thetahatt;
end
if(eta==0)
  tl_etawc=tl_deltasc;
  tl_etawt=tl_deltast;
else
  tl_etawc=tl_eta;
  tl_etawt=tl_eta;
end

% Use stokes 2nd order theory to get vertical fluid velocities
tl_b = -1/r_r2012^2*(1-sqrt(1-r_r2012^2))*tl_r_r2012 ...
       + 1/r_r2012/sqrt(1-r_r2012^2)*r_r2012*tl_r_r2012;
tl_RR = -0.5*tl_b*sin(phi_r2012) ...
        -0.5*b*cos(phi_r2012)*tl_phi_r2012;
tl_worb1c = + pi*tl_Hmo*etawc/(T*h) ...
    + pi*Hmo*tl_etawc/(T*h) ...
    - pi*Hmo*etawc/(T^2*h)*tl_T ...
    - pi*Hmo*etawc/(T*h^2)*tl_h;
tl_worb1t = + pi*tl_Hmo*etawt/(T*h) ...
    + pi*Hmo*tl_etawt/(T*h) ...
    - pi*Hmo*etawt/(T^2*h)*tl_T ...
    - pi*Hmo*etawt/(T*h^2)*tl_h;
tl_worb2c = 2*tl_worb1c*(2*RR-1) ...
    + 2*worb1c*2*tl_RR;
tl_worb2t = 2*tl_worb1t*(2*RR-1) ...
    + 2*worb1t*2*tl_RR;
tl_t1ca = 2*worb1c*tl_worb1c + 32*2*worb2c*tl_worb2c;
tl_t1ta = 2*worb1t*tl_worb1t + 32*2*worb2t*tl_worb2t;
tl_t1cb = tl_worb1c - .5/sqrt(t1ca)*tl_t1ca;
tl_t1tb = tl_worb1t - .5/sqrt(t1ta)*tl_t1ta;
tl_t1c = 2*t1cb/worb2c^2*tl_t1cb ...
         - 2*t1cb^2/worb2c^3*tl_worb2c;
tl_t1t = 2*t1tb/worb2t^2*tl_t1tb ...
         - 2*t1tb^2/worb2t^3*tl_worb2t;
tl_t2cb = -tl_worb1c + 1/sqrt( worb1c^2 + 32*worb2c^2 )*( worb1c*tl_worb1c + 32*worb2c*tl_worb2c );
tl_t2tb = -tl_worb1t + 1/sqrt( worb1t^2 + 32*worb2t^2 )*( worb1t*tl_worb1t + 32*worb2t*tl_worb2t );
tl_t2ca = .125*tl_t2cb/worb2c ...
          - .125*t2cb/worb2c^2*tl_worb2c;
tl_t2ta = .125*tl_t2tb/worb2t ...
          - .125*t2tb/worb2t^2*tl_worb2t;
tl_t2c = -2/sqrt(1-t2ca^2)*tl_t2ca;
tl_t2t = -2/sqrt(1-t2ta^2)*tl_t2ta;
tl_worbc = .125*tl_worb1c*sqrt(t1c) ...
    + .5*.125*worb1c/sqrt(t1c)*tl_t1c ...
    + tl_worb2c*sin(t2c) ...
    + worb2c*cos(t2c)*tl_t2c;
tl_worbt = .125*tl_worb1t*sqrt(t1t) ...
    + .5*.125*worb1t/sqrt(t1t)*tl_t1t ...
    + tl_worb2t*sin(t2t) ...
    + worb2t*cos(t2t)*tl_t2t;

% apply orbital vels to adjust sediment velocity
% wsc=max(ws-worbc,0.001);
% wst=max(ws+worbt,0.001);
tl_wsc = tl_ws - tl_worbc;
tl_wst = tl_ws + tl_worbt;

% sediment load, eqns 23-28.  Note there seems to be a typo in eqns 27-28,
% the bracketed quotient has an issue with dimensions... I changed it into
% what I think is intended.
% Pc = param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc);  % eqn 27
tl_Pc = tl_param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc) ...
        + param.alpha*( ...
            - tl_param.xi*uhatc./c ...
            - param.xi*tl_uhatc./c ...
            + param.xi*uhatc./c^2*tl_c ...
            ).*etawc./(2*(Tc-Tcu)*wsc) ...
        + param.alpha*(1-param.xi*uhatc./c).*tl_etawc./(2*(Tc-Tcu)*wsc) ...
        - param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)^2*( ...
            + 2*(tl_Tc-tl_Tcu)*wsc ...
            + 2*(Tc-Tcu)*tl_wsc );
% Pt = param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst);  % eqn 28
tl_Pt = tl_param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst) ...
        + param.alpha*( ...
            + tl_param.xi*uhatt./c ...
            + param.xi*tl_uhatt./c ...
            - param.xi*uhatt./c^2*tl_c ...
            ).*etawt./(2*(Tt-Ttu)*wst) ...
        + param.alpha*(1+param.xi*uhatt./c).*tl_etawt./(2*(Tt-Ttu)*wst) ...
        - param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)^2*( ...
            + 2*(tl_Tt-tl_Ttu)*wst ...
            + 2*(Tt-Ttu)*tl_wst );
absthetac=abs(thetac);
absthetat=abs(thetat);
tl_absthetac = sign(thetac)*tl_thetac;
tl_absthetat = sign(thetat)*tl_thetat;
if(abs(thetac)>theta_cr)
  tl_Omegac = tl_param.m*(absthetac-theta_cr).^param.n ...
      + param.n*param.m*(absthetac-theta_cr).^(param.n-1)*(tl_absthetac-tl_theta_cr) ...
      + param.m*(absthetac-theta_cr).^param.n*log(absthetac-theta_cr)*tl_param.n;
else
  tl_Omegac=0;
end
if(abs(thetat)>theta_cr)
  tl_Omegat = tl_param.m*(absthetat-theta_cr).^param.n ...
      + param.n*param.m*(absthetat-theta_cr).^(param.n-1)*(tl_absthetat-tl_theta_cr) ...
      + param.m*(absthetat-theta_cr).^param.n*log(absthetat-theta_cr)*tl_param.n;
else
  tl_Omegat=0;
end
if(Pc>1)
  tl_Omegacc = tl_Omegac./Pc ...
      - Omegac./Pc^2*tl_Pc;
else
  tl_Omegacc=tl_Omegac;
end
if(Pt>1)
  tl_Omegatt = tl_Omegat./Pt ...
      - Omegat./Pt^2*tl_Pt;
else
  tl_Omegatt=tl_Omegat;
end
if(Pc<=1)
  tl_Omegact=0;
else
  tl_Omegact = 1./Pc^2*Omegac*tl_Pc ...
      + (1-1./Pc)*tl_Omegac;
end
if(Pt<=1)
  tl_Omegatc=0;
else
  tl_Omegatc = 1./Pt^2*Omegat*tl_Pt ...
      + (1-1./Pt)*tl_Omegat;
end

% transport, eqn 1
% qsc=sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac;
tl_qsc = .5/sqrt(absthetac).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./absthetac.*tl_absthetac ...
         + sqrt(absthetac).*tl_Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./absthetac ...
         + sqrt(absthetac).*Tc.*( ...
             + tl_Omegacc ...
             + tl_Tc./(2*Tcu).*Omegatc ...
             - Tc./(2*Tcu)^2*2*tl_Tcu.*Omegatc ...
             + Tc./(2*Tcu).*tl_Omegatc ...
             ).*thetacx./absthetac ...
         + sqrt(absthetac).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*tl_thetacx./absthetac ...
         - sqrt(absthetac).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./absthetac^2*tl_absthetac;
% qst=sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat;
tl_qst = .5/sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat*tl_absthetat ...
         + sqrt(absthetat)*tl_Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat ...
         + sqrt(absthetat)*Tt*( ...
             + tl_Omegatt ...
             + tl_Tt/(2*Ttu)*Omegact ...
             - Tt/(2*Ttu)^2*2*tl_Ttu*Omegact ...
             + Tt/(2*Ttu)*tl_Omegact ...
             )*thetatx/absthetat ...
         + sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*tl_thetatx/absthetat ...
         - sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat^2*tl_absthetat;
term3 = sqrt((s-1)*g*d50^3)/(1-psed);
tl_term3 = .5/sqrt((s-1)*g*d50^3)*3*(s-1)*g*d50^2*tl_d50/(1-psed);

% Add suspended load contribution, based on energetics model
if(isfield(param,'Cc'))  % OPTION-1

  % qsVdA = (qsc + qst)/T*term3;
  tl_qsVdA = (tl_qsc + tl_qst)/T*term3 ...
        - (qsc + qst)/T^2*term3*tl_T ...
        + (qsc + qst)/T*tl_term3;

  tl_uwmo=tl_uw/1.4;
  % arg_qs2 = (uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)*udelta(1);
  tl_arg_qs2 = 3*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(1/2)*udelta(1).*( ...
      uwmo.*tl_uwmo + udelta(1)*tl_udelta(1) + udelta(2)*tl_udelta(2) ) ...
      + (uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)*tl_udelta(1);
  tl_qs2 = mean(tl_arg_qs2);
  % arg_qs3 = (uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(5/2);
  tl_arg_qs3 = 5*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2).*( ...
      uwmo.*tl_uwmo + udelta(1)*tl_udelta(1) + udelta(2)*tl_udelta(2) );
  tl_qs3 = mean(tl_arg_qs3);
  tl_qsCc = + tl_qs2*eps_s/ws*param.Cc     /(g*(s-1)*(1-psed)) ...
            - qs2*eps_s/ws^2*param.Cc*tl_ws/(g*(s-1)*(1-psed)) ...
            + qs2*eps_s/ws*tl_param.Cc     /(g*(s-1)*(1-psed));
  tl_qsCf = - tl_qs3*eps_s^2/ws^2*param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
            + 2*qs3*eps_s^2/ws^3*tl_ws*param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
            - qs3*eps_s^2/ws^2*tl_param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
            - qs3*eps_s^2/ws^2*param.Cf*tl_tanbeta/(g*(s-1)*(1-psed));
  error('this option should include qbCf, needs to be added')

elseif(~isfield(param,'nosusp') | param.nosusp==0)  % OPTION-2

  if(Omegat>0)
    % Lt=max(eps,fzero(@(L)Omegat*d50*exp(-deltast/L)-0.08*L,.1));
    numst=exp(-deltast/Lt)/0.08;
    denomst=1-Omegat*d50/0.08*deltast/Lt^2*exp(-deltast/Lt);
    dLtdOmegat = numst*d50/denomst;
    dLtdd50 = numst*Omegat/denomst;
    tl_Lt = + dLtdOmegat*tl_Omegat ...
            + dLtdd50*tl_d50;
  else
    tl_Lt=0;
  end
  if(Omegac>0)
    % Lc=max(eps,fzero(@(L)Omegac*d50*exp(-deltasc/L)-0.08*L,.1));
    numsc=exp(-deltasc/Lc)/0.08;
    denomsc=1-Omegac*d50/0.08*deltasc/Lc^2*exp(-deltasc/Lc);
    dLcdOmegac = numsc*d50/denomsc;
    dLcdd50 = numsc*Omegac/denomsc;
    tl_Lc = + dLcdOmegac*tl_Omegac ...
            + dLcdd50*tl_d50;
  else
    tl_Lc=0;
  end
  if(Lc>0)
    % wfracc = (1-exp(-delta/Lc));
    tl_wfracc = - delta/Lc^2*exp(-delta/Lc)*tl_Lc ...
        + 1/Lc*exp(-delta/Lc)*tl_delta;
  else
    tl_wfracc=0;
  end
  if(Lt>0)
    % wfract = (1-exp(-delta/Lt));
    tl_wfract = - delta/Lt^2*exp(-delta/Lt)*tl_Lt ...
        + 1/Lt*exp(-delta/Lt)*tl_delta;
  else
    tl_wfract=0;
  end
  % qsVdA = ( qsc*wfracc + qst*wfract )/T*sqrt((s-1)*g*d50^3)/(1-psed);  % redefined
  tl_qsVdA = (tl_qsc*wfracc + qsc*tl_wfracc + tl_qst*wfract + qst*tl_wfract)/T*term3 ...
      - (qsc*wfracc + qst*wfract)/T^2*term3*tl_T ...
      + (qsc*wfracc + qst*wfract)/T*tl_term3;

  % uwmo=uw/1.4;
  tl_uwmo = + 1/1.4*tl_uw;
  % utot=sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2);  % total current magnitude (time dependent)
  tl_utot = 1./sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).*( ...
      + uwmo     .*tl_uwmo      ...
      + udelta(1)*tl_udelta(1) ...
      + udelta(2)*tl_udelta(2) );
  % Omegai = (1-wfracc)*Omegac*Tc/T + (1-wfract)*Omegat*Tt/T;  % total load above-WBL
  tl_Omegai = ...
      - tl_wfracc*Omegac*Tc/T ...
      + (1-wfracc)*tl_Omegac*Tc/T ...
      + (1-wfracc)*Omegac*tl_Tc/T ...
      - (1-wfracc)*Omegac*Tc/T^2*tl_T ...
      - tl_wfract*Omegat*Tt/T ...
      + (1-wfract)*tl_Omegat*Tt/T ...
      + (1-wfract)*Omegat*tl_Tt/T ...
      - (1-wfract)*Omegat*Tt/T^2*tl_T;
  % Ksdenom1 = + eps_s/ws*mean(utot.^3);
  arg_Ksdenom1 = utot.^3;
  tl_arg_Ksdenom1 = 3*utot.^2.*tl_utot;
  tl_Ksdenom1 = ...
      - eps_s/ws^2*mean(arg_Ksdenom1)*tl_ws ...
      + eps_s/ws*mean(tl_arg_Ksdenom1);
  % Ksdenom2 = - (eps_s/ws)^2*mean(utot.^4)*tanbeta;
  arg_Ksdenom2 = utot.^4;
  tl_arg_Ksdenom2 = 4*utot.^3.*tl_utot;
  tl_Ksdenom2 = + 2*eps_s^2/ws^3*mean(arg_Ksdenom2)*tanbeta*tl_ws ...
      - (eps_s/ws)^2*mean(tl_arg_Ksdenom2)*tanbeta ...
      - (eps_s/ws)^2*mean(arg_Ksdenom2)*tl_tanbeta;
  % Ki = Omegai*(s-1)*g*d50/(Ksdenom1+Ksdenom2);  % chosen s.t. EM-model predicts sediment load = Omegai
  tl_Ki = tl_Omegai*(s-1)*g*d50/(Ksdenom1+Ksdenom2) ...
          + Omegai*(s-1)*g*tl_d50/(Ksdenom1+Ksdenom2) ...
          - Omegai*(s-1)*g*d50/(Ksdenom1+Ksdenom2)^2*(tl_Ksdenom1+tl_Ksdenom2);
  % qsCc = + 1/(g*(s-1)*(1-psed))*eps_s*Ki/ws*mean(arg_qsCc);
  arg_qsCc = utot.^3*udelta(1);
  tl_arg_qsCc = + 3*utot.^2*udelta(1).*tl_utot ...
      + utot.^3*tl_udelta(1);
  tl_qsCc = ...
      + 1/(g*(s-1)*(1-psed))*eps_s*tl_Ki/ws*mean(arg_qsCc) ...
      - 1/(g*(s-1)*(1-psed))*eps_s*Ki/ws^2*mean(arg_qsCc)*tl_ws ...
      + 1/(g*(s-1)*(1-psed))*eps_s*Ki/ws*mean(tl_arg_qsCc);

  % slope-driven transport includes both suspended and bed load.  Use EM-model
  % to predict the ratio between the two
  % Kbdenom1 = + eps_b/tan_phi*mean(utot.^2);
  arg_Kbdenom1 = utot.^2;
  tl_arg_Kbdenom1 = 2*utot.*tl_utot;
  tl_Kbdenom1 = + eps_b/tan_phi*mean(tl_arg_Kbdenom1);
  % Kbdenom2 = - eps_b/tan_phi^2*mean(utot.^2)*tanbeta;
  arg_Kbdenom2 = utot.^2;
  tl_arg_Kbdenom2 = 2*utot.*tl_utot;
  tl_Kbdenom2 = - eps_b/tan_phi^2*mean(tl_arg_Kbdenom2)*tanbeta ...
      - eps_b/tan_phi^2*mean(arg_Kbdenom2)*tl_tanbeta;
  % Omega = Omegac*Tc/T + Omegat*Tt/T;  % total sediment load
  tl_Omega = ...
      + tl_Omegac*Tc/T ...
      + Omegac*tl_Tc/T ...
      - Omegac*Tc/T^2*tl_T ...
      + tl_Omegat*Tt/T ...
      + Omegat*tl_Tt/T ...
      - Omegat*Tt/T^2*tl_T;
  % K = Omega*(s-1)*g*d50/(Kbdenom1+Kbdenom2+Ksdenom1+Ksdenom2);  % chosen s.t. EM-model predicts total sed load = Omega
  tl_K = ...
      + tl_Omega*(s-1)*g*d50/(Kbdenom1+Kbdenom2+Ksdenom1+Ksdenom2) ...
      + Omega*(s-1)*g*tl_d50/(Kbdenom1+Kbdenom2+Ksdenom1+Ksdenom2) ...
      - Omega*(s-1)*g*d50/(Kbdenom1+Kbdenom2+Ksdenom1+Ksdenom2)^2*( tl_Kbdenom1+tl_Kbdenom2+tl_Ksdenom1+tl_Ksdenom2 );
  % qsCf = - slopeFact*K*(eps_s/ws)^2*mean(utot.^5)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven suspended load
  arg_qsCf = utot.^5;
  tl_arg_qsCf = 5*utot.^4.*tl_utot;
  tl_qsCf = ...
      - slopeFact*tl_K*eps_s^2/ws^2*mean(arg_qsCf)*tanbeta/(g*(s-1)*(1-psed)) ...
      + 2*slopeFact*K*eps_s^2/ws^3*mean(arg_qsCf)*tanbeta/(g*(s-1)*(1-psed))*tl_ws ...
      - slopeFact*K*eps_s^2/ws^2*mean(tl_arg_qsCf)*tanbeta/(g*(s-1)*(1-psed)) ...
      - slopeFact*K*eps_s^2/ws^2*mean(arg_qsCf)*tl_tanbeta/(g*(s-1)*(1-psed));
  % qbCf = - slopeFact*K*eps_b/tan_phi^2*mean(utot.^3)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven bed load
  arg_qbCf = utot.^3;
  tl_arg_qbCf = 3*utot.^2.*tl_utot;
  tl_qbCf = ...
      - slopeFact*tl_K*eps_b/tan_phi^2*mean(arg_qbCf)*tanbeta/(g*(s-1)*(1-psed)) ...
      - slopeFact*K*eps_b/tan_phi^2*mean(tl_arg_qbCf)*tanbeta/(g*(s-1)*(1-psed)) ...
      - slopeFact*K*eps_b/tan_phi^2*mean(arg_qbCf)*tl_tanbeta/(g*(s-1)*(1-psed));

else  % OPTION-3, above-WBL transport disabled
  tl_qsVdA = (tl_qsc + tl_qst)/T*term3 ...
      - (qsc + qst)/T^2*term3*tl_T ...
      + (qsc + qst)/T*tl_term3;
  tl_qsCc=0;
  tl_qsCf=0;
  tl_qBCf=0;
end

tl_qs = tl_qsVdA + tl_qsCc + tl_qsCf + tl_qbCf;

% dump all ancillary variables
tl_all.absthetac      =tl_absthetac      ;
tl_all.absthetat      =tl_absthetat      ;
tl_all.ahat           =tl_ahat           ;
tl_all.all            =tl_all            ;
tl_all.alpha          =tl_alpha          ;
tl_all.argc           =tl_argc           ;
tl_all.argc1          =tl_argc1          ;
tl_all.argc2          =tl_argc2          ;
tl_all.argt           =tl_argt           ;
tl_all.argt1          =tl_argt1          ;
tl_all.argt2          =tl_argt2          ;
tl_all.asinarg        =tl_asinarg        ;
tl_all.b              =tl_b              ;
tl_all.c              =tl_c              ;
tl_all.deltasc        =tl_deltasc        ;
tl_all.deltasc        =tl_deltasc        ;
tl_all.deltast        =tl_deltast        ;
tl_all.deltast        =tl_deltast        ;
tl_all.Dstar          =tl_Dstar          ;
tl_all.eta            =tl_eta            ;
tl_all.etawc          =tl_etawc          ;
tl_all.etawt          =tl_etawt          ;
tl_all.fwc            =tl_fwc            ;
tl_all.fwdc           =tl_fwdc           ;
tl_all.fwdt           =tl_fwdt           ;
tl_all.fwt            =tl_fwt            ;
tl_all.Hmo            =tl_Hmo            ;
tl_all.lambda         =tl_lambda         ;
tl_all.meta           =tl_meta           ;
tl_all.mlambda        =tl_mlambda        ;
tl_all.mu             =tl_mu             ;
tl_all.neta           =tl_neta           ;
tl_all.nlambda        =tl_nlambda        ;
tl_all.Omegac         =tl_Omegac         ;
tl_all.Omegacc        =tl_Omegacc        ;
tl_all.Omegact        =tl_Omegact        ;
tl_all.Omegat         =tl_Omegat         ;
tl_all.Omegatc        =tl_Omegatc        ;
tl_all.Omegatt        =tl_Omegatt        ;
tl_all.Pc             =tl_Pc             ;
tl_all.phidc          =tl_phidc          ;
tl_all.phiuc          =tl_phiuc          ;
tl_all.psihat         =tl_psihat         ;
tl_all.psihatc        =tl_psihatc        ;
tl_all.psihatt        =tl_psihatt        ;
tl_all.Pt             =tl_Pt             ;
tl_all.qbCf           =tl_qbCf           ;
tl_all.qs             =tl_qs             ;
tl_all.qsc            =tl_qsc            ;
tl_all.qsCc           =tl_qsCc           ;
tl_all.qsCf           =tl_qsCf           ;
tl_all.qst            =tl_qst            ;
tl_all.qsVdA          =tl_qsVdA          ;
tl_all.RR             =tl_RR             ;
tl_all.streamingEffect=tl_streamingEffect; 
tl_all.T              =tl_T              ;
tl_all.t1c            =tl_t1c            ;
tl_all.t1ca           =tl_t1ca           ;
tl_all.t1cb           =tl_t1cb           ;
tl_all.t1t            =tl_t1t            ;
tl_all.t1ta           =tl_t1ta           ;
tl_all.t1tb           =tl_t1tb           ;
tl_all.t2c            =tl_t2c            ;
tl_all.t2ca           =tl_t2ca           ;
tl_all.t2cb           =tl_t2cb           ;
tl_all.t2t            =tl_t2t            ;
tl_all.t2ta           =tl_t2ta           ;
tl_all.t2tb           =tl_t2tb           ;
tl_all.Tc             =tl_Tc             ;
tl_all.tcr            =tl_tcr            ;
tl_all.Tcu            =tl_Tcu            ;
tl_all.tdc            =tl_tdc            ;
tl_all.term3          =tl_term3          ;
tl_all.thetac         =tl_thetac         ;
tl_all.theta_cr       =tl_theta_cr       ;
tl_all.thetacx        =tl_thetacx        ;
tl_all.thetahatc      =tl_thetahatc      ;
tl_all.thetahatt      =tl_thetahatt      ;
tl_all.thetat         =tl_thetat         ;
tl_all.thetatx        =tl_thetatx        ;
tl_all.Tt             =tl_Tt             ;
tl_all.ttr            =tl_ttr            ;
tl_all.Ttu            =tl_Ttu            ;
tl_all.tuc            =tl_tuc            ;
tl_all.ucrabs         =tl_ucrabs         ;
tl_all.ucrvec         =tl_ucrvec         ;
tl_all.udabs          =tl_udabs          ;
tl_all.uhat           =tl_uhat           ;
tl_all.uhatc          =tl_uhatc          ;
tl_all.uhatt          =tl_uhatt          ;
tl_all.utildecr       =tl_utildecr       ;
tl_all.utildetr       =tl_utildetr       ;
tl_all.utrabs         =tl_utrabs         ;
tl_all.utrvec         =tl_utrvec         ;
tl_all.uw2mean        =tl_uw2mean        ;
tl_all.worb1c         =tl_worb1c         ;
tl_all.worb1t         =tl_worb1t         ;
tl_all.worb2c         =tl_worb2c         ;
tl_all.worb2t         =tl_worb2t         ;
tl_all.worbc          =tl_worbc          ;
tl_all.worbt          =tl_worbt          ;
tl_all.wsc            =tl_wsc            ;
tl_all.wst            =tl_wst            ;
if(isfield(param,'Cc'))  % OPTION-1
  tl_all.uwmo   =tl_uwmo   ;
  tl_all.arg_qs2=tl_arg_qs2; 
  tl_all.qs2    =tl_qs2    ;
  tl_all.arg_qs3=tl_arg_qs3; 
  tl_all.qs3    =tl_qs3    ;
elseif(~isfield(param,'nosusp') | param.nosusp==0)  % OPTION-2
  tl_all.Lt          =tl_Lt          ;
  tl_all.Lc          =tl_Lc          ;
  tl_all.wfracc      =tl_wfracc      ;
  tl_all.wfract      =tl_wfract      ;
  tl_all.uwmo        =tl_uwmo        ;
  tl_all.utot        =tl_utot        ;
  tl_all.Omegai      =tl_Omegai      ;
  tl_all.Omega       =tl_Omega       ;
  tl_all.arg_Ksdenom1=tl_arg_Ksdenom1; 
  tl_all.arg_Ksdenom2=tl_arg_Ksdenom2; 
  tl_all.Ksdenom1    =tl_Ksdenom1    ;
  tl_all.Ksdenom2    =tl_Ksdenom2    ;
  tl_all.Ki          =tl_Ki          ;
  tl_all.arg_Kbdenom1=tl_arg_Kbdenom1; 
  tl_all.arg_Kbdenom2=tl_arg_Kbdenom2; 
  tl_all.Kbdenom1    =tl_Kbdenom1    ;
  tl_all.Kbdenom2    =tl_Kbdenom2    ;
  tl_all.K           =tl_K           ;
  tl_all.arg_qsCc    =tl_arg_qsCc    ;
  tl_all.arg_qsCf    =tl_arg_qsCf    ;
  tl_all.arg_qbCf    =tl_arg_qbCf    ;
end
if(~isfield(param,'streamingType') || param.streamingType=='v')
  tl_all.fwd=tl_fwd;
elseif(param.streamingType=='n')
  tl_all.f25    =tl_f25    ;
  tl_all.theta25=tl_theta25;
  tl_all.r      =tl_r      ;
  tl_all.fws    =tl_fws    ;
end
tl_all.tauwRe=tl_tauwRe;

% TEST-CODE: override output variable
% eval(['tl_qs = tl_' outvar ';']);

end  % end of main function
