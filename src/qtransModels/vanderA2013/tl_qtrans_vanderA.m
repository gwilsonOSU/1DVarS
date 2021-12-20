function tl_qs=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd,outvar)
%
% TL code for qtrans_vanderA.m
%

% this wrapper loop serves to handle vector inputs
nx=length(tl_h);
tl_Q=zeros(nx,1);
for i=1:nx
  tl_qs(i)= ...
      tl_qtrans_vanderA_main(tl_d50(i),tl_d90(i),tl_h(i),tl_Hrms(i),tl_kabs(i),...
                             tl_omega,tl_udelta(i,:),tl_ws(i),tl_Aw(i),tl_Sw(i),tl_Uw(i),tl_param,bkgd(i),outvar);
end
tl_qs=tl_qs(:);

end  % end of wrapper function, start of main function

function tl_qs=tl_qtrans_vanderA_main(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd,outvar)

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
delta       =bkgd.delta       ;
deltasc     =bkgd.deltasc     ;
deltast     =bkgd.deltast     ;
eta         =bkgd.eta         ;
etawc       =bkgd.etawc       ;
etawt       =bkgd.etawt       ;
fd          =bkgd.fd          ;
fw          =bkgd.fw          ;
fwc         =bkgd.fwc         ;
fwd         =bkgd.fwd         ;
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
worbc       =bkgd.worbc       ;
worbt       =bkgd.worbt       ;
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

% derived params
tl_Hmo=tl_Hrms*1.4;
tl_c=-omega/kabs^2*tl_kabs + tl_omega/kabs;

% intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% relevant to van Der A
[tl_uw,tl_r_r2012,tl_phi_r2012]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,omega*t,uwave_wksp);
uw2mean=mean(uw.^2);
tl_uw2mean=mean(2*uw.*tl_uw);
tl_uhat=.5./sqrt(2*uw2mean)*2.*tl_uw2mean;

% timing of wave velocity direction change, crest, and trough, based on
% Ruessink et al 2012.
asinarg=-r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2));
tl_asinarg = -tl_r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)) ...
    - r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))*tl_phi_r2012 ...
    + r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)).^2.*( ...
        .5./sqrt(1-r_r2012^2)*(-2*r_r2012*tl_r_r2012) );

% Tc=asin(asinarg)/omega;
tl_Tc = 1/sqrt(1-asinarg^2)*tl_asinarg/omega ...
        - asin(asinarg)/omega^2*tl_omega;
tl_Tt=-tl_Tc;
tl_T=tl_Tc+tl_Tt;
tl_Tcu = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,...
                                      omega,r_r2012,phi_r2012,Tcu);
tl_Ttu = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,...
                                      omega,r_r2012,phi_r2012,Ttu);
tl_Ttu=tl_Tt-tl_Tc;  % duration of deceleration under trough

% for crest velocities, best I can do without being too fancy and screwing
% up the AD code is just to consider the TL velocity at the location of the
% NL model crest/trough
[~,ic]=max(uw);
[~,it]=min(uw);
tl_uhatc = tl_uw(ic);
tl_uhatt = -tl_uw(it);

% critical shields param, Soulsby
tl_Dstar=(g*(s-1)/nu^2)^(1/3)*tl_d50;
% theta_cr=.3/(1+1.2*Dstar)+0.055*(1-exp(-.02*Dstar));
tl_theta_cr = .3/(1+1.2*Dstar)^2*(1.2*tl_Dstar) ...
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
    tl_vanderA_shields(tl_d50,tl_d90,tl_udabs,tl_uhat,...
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
tl_fwd = tl_alpha*fd ...
         + alpha*tl_fd ...
         - tl_alpha*fw ...
         + (1-alpha)*tl_fw;
tl_tauwRe = .5*alphaw*( ...
    tl_fwd*uhat^3/c ...
    + 3*fwd*uhat^2/c*tl_uhat ...
    - fwd*uhat^3/c^2*tl_c );
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

% % Use stokes 2nd order theory to get vertical fluid velocities and correct
% % settling velocity.  Follows Malarkey & Davies (2012).  Some of this is
% % simply ported from COAWST code.  NOTE, it is a bit odd to do hydrodynamic
% % calculations in qtrans_vanderA (by convention, hydro stuff is ideally done
% % by callers of this function), but the code can't be easily moved outside
% % because it uses ripple height as input.
% % b=1/r_r2012*(1-sqrt(1-r_r2012^2));  % see MD2012, line after eqn 13b
% tl_b = - 1/r_r2012^2*(1-sqrt(1-r_r2012^2))*tl_r_r2012 ...
%     + 1/r_r2012/sqrt(1-r_r2012^2)*r_r2012*tl_r_r2012;
% % RR=0.5*(1+b*sin(-phi_r2012));  % see MD2012, eqn 17, note their phi convention is negated
% tl_RR = 0.5*sin(-phi_r2012)*tl_b ...
%         - 0.5*b*cos(-phi_r2012)*tl_phi_r2012;
% % worb1c=pi*Hmo*etawc/(T*h);
% tl_worb1c = pi*tl_Hmo*etawc/(T*h) ...
%     + pi*Hmo*tl_etawc/(T*h) ...
%     - pi*Hmo*etawc/(T^2*h)*tl_T ...
%     - pi*Hmo*etawc/(T*h^2)*tl_h;
% % worb1t=pi*Hmo*etawt/(T*h);
% tl_worb1t = pi*tl_Hmo*etawt/(T*h) ...
%     + pi*Hmo*tl_etawt/(T*h) ...
%     - pi*Hmo*etawt/(T^2*h)*tl_T ...
%     - pi*Hmo*etawt/(T*h^2)*tl_h;
% % worb2c=2*worb1c*(2*RR-1);
% tl_worb2c = 2*tl_worb1c*(2*RR-1) ...
%     + 2*worb1c*2*tl_RR;
% % worb2t=2*worb1t*(2*RR-1);
% tl_worb2t = 2*tl_worb1t*(2*RR-1) ...
%     + 2*worb1t*2*tl_RR;
% tl_t1ca = 2*worb1c*tl_worb1c + 2*32*worb2c*tl_worb2c;
% tl_t1ta = 2*worb1t*tl_worb1t + 2*32*worb2t*tl_worb2t;
% tl_t1cb = tl_worb1c - .5/sqrt(t1ca)*tl_t1ca;
% tl_t1tb = tl_worb1t - .5/sqrt(t1ta)*tl_t1ta;
% tl_t1c = 2*t1cb/worb2c^2*tl_t1cb - 2*t1cb^2/worb2c^3*tl_worb2c;
% tl_t1t = 2*t1tb/worb2t^2*tl_t1tb - 2*t1tb^2/worb2t^3*tl_worb2t;
% tl_t2cb = -tl_worb1c + 1/sqrt( worb1c^2 + 32*worb2c^2 )*( worb1c*tl_worb1c + 32*worb2c*tl_worb2c );
% tl_t2tb = -tl_worb1t + 1/sqrt( worb1t^2 + 32*worb2t^2 )*( worb1t*tl_worb1t + 32*worb2t*tl_worb2t );
% tl_t2ca = .125*tl_t2cb/worb2c - .125*t2cb/worb2c^2*tl_worb2c;
% tl_t2ta = .125*tl_t2tb/worb2t - .125*t2tb/worb2t^2*tl_worb2t;
% tl_t2c = 2/sqrt(1-t2ca^2)*tl_t2ca;
% tl_t2t = 2/sqrt(1-t2ta^2)*tl_t2ta;
% tl_worbc = .125*tl_worb1c*sqrt(t1c) ...
%     + .5*.125*worb1c/sqrt(t1c)*tl_t1c ...
%     + tl_worb2c*sin(t2c) ...
%     + worb2c*cos(t2c)*tl_t2c;
% tl_worbt = .125*tl_worb1t*sqrt(t1t) ...
%     + .5*.125*worb1t/sqrt(t1t)*tl_t1t ...
%     + tl_worb2t*sin(t2t) ...
%     + worb2t*cos(t2t)*tl_t2t;

% TODO: for now, don't use the above block until I get the TL-AD symmetry fixed
tl_worbc=0;
tl_worbt=0;

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
qsc=sqrt(absthetac).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./absthetac;
qst=sqrt(absthetat).*Tt.*(Omegatt+Tt./(2*Ttu).*Omegact).*thetatx./absthetat;
term3=sqrt((s-1)*g*d50^3)/(1-psed);
qs = ( qsc + qst )./T.*term3;

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
% qs = (qsc + qst)/T*term3;
tl_qs = (tl_qsc + tl_qst)/T*term3 ...
        - (qsc + qst)/T^2*term3*tl_T ...
        + (qsc + qst)/T*tl_term3;

% TEST-CODE: override output variable
eval(['tl_qs = tl_' outvar ';']);

end  % end of main function
