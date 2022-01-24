function [qs,workspc,Omegatc]=qtrans_vanderA_onewave(uw,d50,d90,wsc,wst,udelta,delta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,~,~,tanbeta,param,Omegatc)
%
% [qs,workspc]=qtrans_vanderA(uw,d50,d90,wsc,wst,udelta,delta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param,Omegatc)
%
% Calculates transport following van der A (2013).  This version is for a
% single wave with measured wave shape parameters (uhat,uhatc,T,Tc,etc etc).
% If you only have bulk wave parameters (Hrms,omega), you can instead use
% qtrans_vanderA.m.
%
% INPUTS:
%
% uw : time series of cross-shore component of velocity for complete wave
%      period.  Is used for estimation of suspended sediment above the WBL.
% d50   : in meters
% d90   : in meters 
% wsc,wst : particle settling velocities under crest and trough, m/s,
%           potentially including a correction for orbital velocity.  Note:
%           If you don't have measurements these can be estimated from wave
%           parameters instead; refer to the implementation in
%           qtrans_vanderA.m as an example.
% udelta: near-bed mean velocity, 2x1 vector, m/s.  Note: wave angle is
%         assumed to be zero in this code.  To incorporate a nonzero wave
%         angle, you should rotate udelta to align with the wave direction.
%         An example of doing this can be found in hydroSedModel.m provided
%         with the 1DVarS code repository.
% uhat=sqrt(2*mean(uw.^2)) : rms wave velocity for full wave cycle, eqn 8, m/s
% uhatc=max(uw) : max velocity under crest, see text below eqn 7, m/s
% uhatt=min(uw) : max velocity under trough (positive valued), m/s
% T   : full wave period, s
% Tc  : duration of crest, s
% Tt  : duration of trough, s
% Ttu : duration of deceleration under trough, s
% Tcu : duration of acceleration under crest, s
% c   : wave phase speed, m/s (Note: Set param.xi for flow-tunnel experiments, in which case the value of c is irrelevant to the result)
% eta,lambda: Ripple height and wavelength.  Set eta==0 to neglect ripple effects.
% tanbeta : tangent of beach slope.  Used for calculating slope-driven suspended transport above the WBL.
% param.{alpha,xi,m,n}: tuning parameters struct
%        alpha: phase lag effect, eqns (24)-(28).  Default 8.2
%        xi   : phase lag effect, eqns (24)-(28).  Default 1.7.  O(1) according to Kranenburg (2013)
%        m    : MPM leading coefficient, eqn (36).  Default 11.0
%        n    : MPM exponent, eqn (36).  Default 1.2
%
% OUTPUTS:
%
% qs: total sediment transport flux in m2/s.  This code does not include any
% correction for bed porosity, so qs is volume of sediment per unit time per
% unit alongshore distance.
%
% workspc: this is intended for passing internally to TL-AD codes, not for
% end-user.
%
%     If you wish to get user access to output variables other than q, it
%     would probably be better to add them to the list of outputs rather
%     than attempt to extract them from the 'workspc' struct.
%
% !!!!IMPORTANT NOTE, TODO!!!!
%
% This code is mostly copied from qtrans_vanderA.m, but this version takes
% wave shape parameters (uhat,uhatc,uhatt,T,Tc,Tt,Tcu,Ttu) as inputs, while
% qtrans_vanderA.m calculates wave shape from a model based on bulk wave
% parameters (Hrms,omega).  Ideally, wave shape parameters would be the
% default way of running the model, and when needed the wave shape could be
% calculated from bulk wave parameters externally as a separate function.
% The ONLY reason for not doing so is the calculation of differential
% settling velocities (wsc,wst) from bulk parameters is tied up with the vDA
% model's ripple calculation, making code separation difficult.  In the
% future it would be good to refactor somehow; the current situation of
% having mostly-duplicate codes for handling two different input scenarios
% is not good.
%

physicalConstants;

% critical shields param, Soulsby
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
psihatc=(1.27*uhatc)^2/((s-1)*g*d50);  % see Appendix B.  factor 1.27 is for irregular flow
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

% shields parameter related parameters.  Requires solving a 5-eqn nonlinear
% system, so this is done in its own code.
udabs=sqrt(udelta(1)^2+udelta(2)^2);
[theta_av,ksd,ksw,fd,fw,...
          branch_A1,branch_A4,branch_A5] = ...
    vanderA_shields(d50,d90,udabs,uhat,delta,mu,eta,lambda,ahat);

% other BBL derived parameters
alpha=udabs./(udabs+uhat);  % eqn 19
c1=2.6;  % cited in text as van der A et al. (2011)
if(ahat/ksw>1.587)  % eqn 21
  fwc=.00251*exp(5.21*((2*Tcu/Tc)^c1*ahat/ksw)^-.19);
  fwt=.00251*exp(5.21*((2*Ttu/Tt)^c1*ahat/ksw)^-.19);
else
  fwc=0.3;
  fwt=0.3;
end
fwdc=alpha*fd+(1-alpha)*fwc;  % eqn 18
fwdt=alpha*fd+(1-alpha)*fwt;  % eqn 18
ucrvec=utildecr*[+1 0]+udelta;  % eqn 12.  Assumes x is wave direction.
utrvec=utildetr*[-1 0]+udelta;  % eqn 13
ucrabs=sqrt(ucrvec(1)^2+ucrvec(2)^2);
utrabs=sqrt(utrvec(1)^2+utrvec(2)^2);
thetac=.5*fwdc.*ucrabs.^2/((s-1)*g*d50);  % eqn 17
thetat=.5*fwdt.*utrabs.^2/((s-1)*g*d50);  % eqn 17

% boundary layer streaming, either with VDA13 model or Nielsen
alphaw = 4/(3*pi);
if(~isfield(param,'streamingType') || param.streamingType=='v')
  f25    =nan;  % not used
  theta25=nan;  % not used
  r      =nan;  % not used
  fws    =nan;  % not used
  fwd = alpha*fd+(1-alpha)*fw;
  tauwRe = fwd*alphaw*uhat^3/2/c;
elseif(param.streamingType=='n')
  fwd=nan;  % not used
  f25 = exp(5.5*(2.5*d50/ahat)^.2-6.3);
  theta25 = 0.5*f25*(ahat*omega)^2/((s-1)*g*d50);
  r = 170*sqrt(max(0,theta25-0.05))*d50;
  if(lambda>0)
    r = r + 8*eta^2/lambda;
  end
  fws = exp(5.5*(r/ahat)^.2-6.3);
  tauwRe = fws*alphaw*uhat^3/2/c;
  % tauRe = 1/T/2*fws*trapz(t,abs(uw).^3)/c;  % same but integrate uw(t)
elseif(param.streamingType=='0')
  f25    =nan;  % not used
  theta25=nan;  % not used
  r      =nan;  % not used
  fws    =nan;  % not used
  fwd    =nan;  % not used
  tauwRe=0;
else
  error('must provide param.streamingType as either ''n'', ''v'', or ''0''')
end
streamingEffect = tauwRe/((s-1)*g*d50);  % eqns 15 and 22

% apply streaming to bottom stress
thetacx = abs(thetac)*ucrvec(1)/ucrabs + streamingEffect;  % eqn 15
thetatx = abs(thetat)*utrvec(1)/utrabs + streamingEffect;  % eqn 15

% sheet flow layer thickness, Appendix C
thetahatc=.5*fwdc.*uhatc.^2/((s-1)*g*d50);  % eqn C2
thetahatt=.5*fwdt.*uhatt.^2/((s-1)*g*d50);  % eqn C2
if(d50<=.15e-3)  % eqn C1
  deltasc=(d50)*25*thetahatc;
  deltast=(d50)*25*thetahatt;
elseif(.15e-3<d50&d50<.2e-3)
  deltasc=(d50)*(25-12*(d50-.15e-3)/.05e-3);
  deltast=(d50)*(25-12*(d50-.15e-3)/.05e-3);
else
  deltasc=(d50)*13*thetahatc;
  deltast=(d50)*13*thetahatt;
end
if(eta==0)
  etawc=deltasc;
  etawt=deltast;
else
  etawc=eta;
  etawt=eta;
end

% sediment load, eqns 23-28.  Note there seems to be a typo in eqns 27-28,
% the bracketed quotient has an issue with dimensions... I changed it into
% what I think is intended.
Pc = param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc);  % eqn 27
Pt = param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst);  % eqn 28
if(abs(thetac)>theta_cr)
  Omegac=param.m*(abs(thetac)-theta_cr).^param.n;  % eqn 2
else
  Omegac=0;
end
if(abs(thetat)>theta_cr)
  Omegat=param.m*(abs(thetat)-theta_cr).^param.n;  % eqn 2
else
  Omegat=0;
end

if(Pc>1)
  Omegacc=Omegac./Pc; % eqn 23
else
  Omegacc=Omegac;
end
if(Pt>1)
  Omegatt=Omegat./Pt;  % eqn 25
else
  Omegatt=Omegat;
end
if(Pc<=1)
  Omegact=0;
else
  Omegact=(1-1./Pc)*Omegac;  % eqn 24
end
if exist('Omegatc','var') == 0 % if Omegatc does not exist as an input, calculate Omegatc according to equation 26
    if(Pt<=1)
      Omegatc=0;
    else
      Omegatc=(1-1./Pt)*Omegat;  % eqn 26
    end
else % otherwise, use existing value of Omegatc (presumably from previous wave)
end

% transport, eqn 1
qsc=sqrt(abs(thetac)).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./abs(thetac);
qst=sqrt(abs(thetat)).*Tt.*(Omegatt+Tt./(2*Ttu).*Omegact).*thetatx./abs(thetat);

% Add suspended-load contribution from currents and slope-effect, based on
% energetics model.  This represents transport by sediments suspended above
% the wave boundary layer (WBL) which are not included in the original
% model.
%
% OPTION-1: If 'Cc' is included in params struct, then use explicit
% calculation with EM model.  Code borrowed from qtrans_dubarbier.m.  Note
% the dubarbier code is tuned for rms wave velocity, while VDA uses
% significant wave velocity, so I have to divide uw by 1.4.  The
% "wave-driven" part of the suspended transport is omitted since arguably
% VDA already includes its contribution, and its value is very similar to
% the VDA prediction.
%
% OPTION-2: If 'Cc' is not in the params struct, then the suspended-load
% contribution is based on boundary layer height calculated by VDA model.
% See further info in comments below.
%
if(isfield(param,'Cc'))  % OPTION-1

  eps_s=0.015;
  uwmo=uw/1.4;
  ws=.5*(wsc+wst);
  qs2=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^3*udelta(1));
  qs3=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^5);
  qsCc = + qs2*eps_s/ws*param.Cc                 /(g*(s-1)*(1-psed));
  qsCf = - qs3*eps_s/ws*param.Cf*eps_s*tanbeta/ws/(g*(s-1)*(1-psed));
  qsVdA = ( qsc + qst )./T.*sqrt((s-1)*g*d50^3)/(1-psed);

else  % OPTION-2

  % calculate boundary layer e-folding length based on definition of delta_s
  % by Dohmen- Janssen (1999), c(delta_s) = 0.08.  Then use the calculated
  % e-folding length to determine the fraction of total sediment load
  % (Omegac + Omegat) occurring above the WBL, z>delta.  The above-WBL
  % fraction will be removed from the wave-driven transport calculated
  % previously, and instead used for suspended transport above WBL.
  Lt=max(eps,fzero(@(L)Omegat*d50*exp(-deltast/L)-0.08*L,.1));
  Lc=max(eps,fzero(@(L)Omegac*d50*exp(-deltasc/L)-0.08*L,.1));
  wfracc = (1-exp(-delta/Lc));
  wfract = (1-exp(-delta/Lt));
  qsVdA = ( qsc*wfracc + qst*wfract )./T.*sqrt((s-1)*g*d50^3)/(1-psed);

  % current- and slope-driven transport are partitioned according to
  % energetics model (i.e., "wave stirring" doing work against gravitational
  % settling).  The two effects are assumed share a common drag coefficient,
  % (K = Cc (==Cf) in Dubarbier's notation).
  eps_s=0.015;
  uwmo=uw/1.4;
  ws=.5*(wsc+wst);
  Omegas = (1-wfracc)*Omegac*Tc/T + (1-wfract)*Omegat*Tt/T;  % total load above-WBL
  utot=sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2);  % total current magnitude (time dependent)
  Kdenom1 = + eps_s/ws*mean(utot.^3);
  Kdenom2 = - (eps_s/ws)^2*mean(utot.^4)*tanbeta;
  Kdenom=Kdenom1+Kdenom2;
  K = Omegas*(s-1)*g*d50/(Kdenom1+Kdenom2);  % chosen s.t. EM-model predicts sediment load = Omegas
  qsCc = + K*eps_s/ws*mean(utot.^3*udelta(1))/(g*(s-1)*(1-psed));    % current-driven
  qsCf = - K*(eps_s/ws)^2*mean(utot.^5)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven

end
qs = qsVdA + qsCc + qsCf;

% save all relevant variables for passing to TL model.  Note this runs
% faster if the variables are hard-coded, don't use eval
if(nargout>1)
  workspc=struct;
  workspc.Dstar        =Dstar          ;
  % workspc.Hmo          =Hmo            ;
  % workspc.Hrms         =Hrms           ;
  workspc.Omegac       =Omegac         ;
  workspc.Omegacc      =Omegacc        ;
  workspc.Omegact      =Omegact        ;
  workspc.Omegat       =Omegat         ;
  workspc.Omegatc      =Omegatc        ;
  workspc.Omegatt      =Omegatt        ;
  workspc.Pc           =Pc             ;
  workspc.Pt           =Pt             ;
  workspc.T            =T              ;
  workspc.Tc           =Tc             ;
  workspc.Tcu          =Tcu            ;
  workspc.Tt           =Tt             ;
  workspc.Ttu          =Ttu            ;
  workspc.ahat         =ahat           ;
  workspc.alpha        =alpha          ;
  workspc.alphaw       =alphaw         ;
  workspc.branch_A1    =branch_A1      ;
  workspc.branch_A4    =branch_A4      ;
  workspc.branch_A5    =branch_A5      ;
  workspc.c            =c              ;
  workspc.d50          =d50            ;
  workspc.d90          =d90            ;
  workspc.deltasc      =deltasc        ;
  workspc.deltast      =deltast        ;
  workspc.eta          =eta            ;
  workspc.etawc        =etawc          ;
  workspc.etawt        =etawt          ;
  workspc.fd           =fd             ;
  workspc.fw           =fw             ;
  workspc.fwc          =fwc            ;
  workspc.fwd          =fwd            ;
  workspc.fwdc         =fwdc           ;
  workspc.fwdt         =fwdt           ;
  workspc.fwt          =fwt            ;
  workspc.tauwRe       =tauwRe         ;
  % workspc.h            =h              ;
  % workspc.kabs         =kabs           ;
  workspc.ksd          =ksd            ;
  workspc.ksw          =ksw            ;
  workspc.lambda       =lambda         ;
  workspc.meta         =meta           ;
  workspc.mlambda      =mlambda        ;
  workspc.mu           =mu             ;
  workspc.neta         =neta           ;
  workspc.nlambda      =nlambda        ;
  % workspc.nt           =nt             ;
  % workspc.omega        =omega          ;
  workspc.param        =param          ;
  % workspc.phi_r2012    =phi_r2012      ;
  % workspc.r_r2012      =r_r2012        ;
  workspc.psed         =psed           ;
  workspc.psihat       =psihat         ;
  workspc.psihatc      =psihatc        ;
  workspc.psihatt      =psihatt        ;
  workspc.qs           =qs             ;
  workspc.qsc          =qsc            ;
  workspc.qst          =qst            ;
  % workspc.t            =t              ;
  workspc.theta_av     =theta_av       ;
  workspc.theta_cr     =theta_cr       ;
  workspc.thetac       =thetac         ;
  workspc.thetacx      =thetacx        ;
  workspc.thetatx      =thetatx        ;
  workspc.thetahatc    =thetahatc      ;
  workspc.thetahatt    =thetahatt      ;
  workspc.thetat       =thetat         ;
  workspc.ucrabs       =ucrabs         ;
  workspc.ucrvec       =ucrvec         ;
  workspc.udabs        =udabs          ;
  workspc.udelta       =udelta         ;
  workspc.delta        =delta          ;
  workspc.uhat         =uhat           ;
  workspc.uhatc        =uhatc          ;
  workspc.uhatt        =uhatt          ;
  workspc.utildecr     =utildecr       ;
  workspc.utildetr     =utildetr       ;
  workspc.utrabs       =utrabs         ;
  workspc.utrvec       =utrvec         ;
  workspc.uw           =uw             ;
  workspc.ws           =ws             ;
  workspc.wsc          =wsc            ;
  workspc.wst          =wst            ;
  % workspc.worbc        =worbc          ;
  % workspc.worbt        =worbt          ;
  % workspc.uwave_wksp   =uwave_wksp     ;
  % workspc.Aw           =Aw             ;
  % workspc.Sw           =Sw             ;
  % workspc.Uw           =Uw             ;
  workspc.c1           =c1             ;
  % workspc.b            =b              ;
%   workspc.RR           =RR             ;
%   workspc.worb1c       =worb1c         ;
%   workspc.worb1t       =worb1t         ;
%   workspc.worb2c       =worb2c         ;
%   workspc.worb2t       =worb2t         ;
%   workspc.t1ca         =t1ca           ;
%   workspc.t1ta         =t1ta           ;
%   workspc.t1cb         =t1cb           ;
%   workspc.t1tb         =t1tb           ;
%   workspc.t1c          =t1c            ;
%   workspc.t1t          =t1t            ;
%   workspc.t2cb         =t2cb           ;
%   workspc.t2tb         =t2tb           ;
%   workspc.t2ca         =t2ca           ;
%   workspc.t2ta         =t2ta           ;
%   workspc.t2c          =t2c            ;
%   workspc.t2t          =t2t            ;
  % workspc.worbc        =worbc          ;
  % workspc.worbt        =worbt          ;
  workspc.f25          =f25            ;
  workspc.theta25      =theta25        ;
  workspc.r            =r              ;
  workspc.fws          =fws            ;
%   workspc.phiuc        =phiuc          ;
%   workspc.phidc        =phidc          ;
%   workspc.icu_guess    =icu_guess      ;
%   workspc.itu_guess    =itu_guess      ;

  % accounting of suspended seds above WBL
  workspc.eps_s        =eps_s          ;
  workspc.uwmo         =uwmo           ;
  workspc.tanbeta      =tanbeta        ;
  if(isfield(param,'Cc'))
    workspc.qs2          =qs2            ;
    workspc.qs3          =qs3            ;
  else
    workspc.Lt     =Lt      ;
    workspc.Lc     =Lc      ;
    workspc.wfracc =wfracc  ;
    workspc.wfract =wfract  ;
    workspc.utot   =utot    ;
    workspc.Omegas =Omegas  ;
    workspc.Kdenom1=Kdenom1 ;
    workspc.Kdenom2=Kdenom2 ;
    workspc.Kdenom =Kdenom  ;
    workspc.K      =K       ;
  end
  workspc.qsVdA        =qsVdA          ;
  workspc.qsCc         =qsCc           ;
  workspc.qsCf         =qsCf           ;

end
