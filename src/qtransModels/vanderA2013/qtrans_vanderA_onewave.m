function [qs,workspc,Omegatc]=qtrans_vanderA_onewave(uw,uhat,uhatc,uhatt,T,Tc,Tt,Tcu,Ttu,c,d50,d90,h,tanbeta,udelta,delta,wsc,wst,param,Omegatc)%,outvar)
%
%
% Calculates transport following van der A (2013).  This function is the
% same as qtrans_vanderA.m, except instead of requiring bulk wave parameters
% it takes individual wave parameters.  (If you were to enter individual
% parameters that happen to coincide with those of the bulk wave parameters,
% you'd get the same result as qtrans_vanderA.m.)
%
% INPUTS:
%
% All inputs are as described in qtrans_vanderA.m (note the order of
% paramters is different), WITH THE EXCEPTTION OF...
%
% uw = wave velocity time series, for full wave phase, m/s
% uhat = sqrt(2*mean(uw.^2)) : rms wave velocity for full wave cycle, eqn 8, m/s
% uhatc = max(uw) : max velocity under crest, see text below eqn 7, m/s
% uhatt = -min(uw) : max velocity under trough (positive valued), m/s
% T   = full wave period, s
% Tc  = duration of crest, s
% Tt  = duration of trough, s
% Tcu = duration of acceleration under trough, s
% Ttu = duration of deceleration under crest, s
% c = wave phase speed, m/s (Note: If simulating a flow-tunnel experiment, c
%                            is irrelevant.  In that case you can set
%                            param.xi=0 and c=inf.)
% wsc = settling velocity under wave crest, m/s
% wst = settling velocity under wave trough, m/s
%
% OUTPUTS:
%
% All outputs are as described in qtrans_vanderA.m
%

if(isfield(param,'Cc'))
  warning(['Explicit EEM-based suspended-load transport is being depreciated at time of writing.  ' ...
           'Instead, if you remove Cc from params struct, the suspended-load contribution ' ...
           'will be calculated automatically']);
end

% this wrapper loop serves to handle vector inputs
nx=length(h);
for i=1:nx
  if exist('Omegatc') == 0
    [qs(i),workspc(i),Omegatc(i)]=qtrans_vanderA_onewave_main(uw(:,i),uhat(i),uhatc(i),uhatt(i),T(i),Tc(i),Tt(i),Tcu(i),Ttu(i),c(i),d50(i),d90(i),h(i),tanbeta(i),udelta(i,:),delta(i),wsc(i),wst(i),param(i));%,outvar)
  else
    [qs(i),workspc(i),Omegatc(i)]=qtrans_vanderA_onewave_main(uw(:,i),uhat(i),uhatc(i),uhatt(i),T(i),Tc(i),Tt(i),Tcu(i),Ttu(i),c(i),d50(i),d90(i),h(i),tanbeta(i),udelta(i,:),delta(i),wsc(i),wst(i),param(i),Omegatc(i));%,outvar)
  end
end
qs=qs(:);
workspc=workspc(:);

end  % end of wrapper function, start of main function

function [qs,workspc,localOmegatc]=qtrans_vanderA_onewave_main(uw,uhat,uhatc,uhatt,T,Tc,Tt,Tcu,Ttu,c,d50,d90,h,tanbeta,udelta,delta,wsc,wst,param,Omegatc)%,outvar)

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
eta=0;  % TEST: disable ripples

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
if(Pt<=1)
  localOmegatc=0;
else
  localOmegatc=(1-1./Pt)*Omegat;  % eqn 26
end

% handle optional passing of Omegatc from one wave to the next: If Omegatc
% is provided as input, use it here.  If not provided as input, assume waves
% are periodic hence the previous wave has Omegatc==localOmegatc.  Note,
% regardless the localOmegatc will always be provided as output.
if(~exist('Omegatc'))
  Omegatc=localOmegatc;
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
  qs2=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^3*udelta(1));
  qs3=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^5);
  ws=.5*(wsc+wst);
  qsCc = + qs2*eps_s/ws*param.Cc                 /(g*(s-1)*(1-psed));
  qsCf = - qs3*eps_s/ws*param.Cf*eps_s*tanbeta/ws/(g*(s-1)*(1-psed));
  qsVdA = ( qsc + qst )./T.*sqrt((s-1)*g*d50^3)/(1-psed);

elseif(~isfield(param,'nosusp') || param.nosusp==0)  % OPTION-2

  % calculate boundary layer e-folding length based on definition of delta_s
  % by Dohmen- Janssen (1999), c(delta_s) = 0.08.  Then use the calculated
  % e-folding length to determine the fraction of total sediment load
  % (Omegac + Omegat) occurring above the WBL, z>delta.  The above-WBL
  % fraction will be removed from the wave-driven transport calculated
  % previously, and instead used for suspended transport above WBL.
  if(Omegat>0)
    Lt=max(eps,fzero(@(L)Omegat*d50*exp(-deltast/L)-0.08*L,.1));
  else
    Lt=0;
  end
  if(Omegac>0)
    Lc=max(eps,fzero(@(L)Omegac*d50*exp(-deltasc/L)-0.08*L,.1));
  else
    Lc=0;
  end
  if(Lc>0)
    wfracc = (1-exp(-delta/Lc));  % above-WBL fraction
  else
    wfracc=0;
  end
  if(Lt>0)
    wfract = (1-exp(-delta/Lt));  % above-WBL fraction
  else
    wfract=0;
  end
  qsVdA = ( qsc*wfracc + qst*wfract )./T.*sqrt((s-1)*g*d50^3)/(1-psed);  % total transport within WBL

  % EM-model variables
  eps_s=0.015;
  eps_b=0.135;
  tan_phi=0.63;  % sediment internal friction angle
  uwmo=uw/1.4;
  utot=sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2);  % total current magnitude (time dependent)
  slopeFact=0.25;  % reduction factor for slope-driven transport efficiency

  % above-WBL transport is according to energetics model (i.e., "wave
  % stirring" doing work against gravitational settling)
  Omegai = (1-wfracc)*Omegac*Tc/T + (1-wfract)*Omegat*Tt/T;  % total load above-WBL
  ws=.5*(wsc+wst);
  Ksdenom1 = + eps_s/ws*mean(utot.^3);
  Ksdenom2 = - slopeFact*(eps_s/ws)^2*mean(utot.^4)*tanbeta;
  Ki = Omegai*(s-1)*g*d50/(Ksdenom1+Ksdenom2);  % chosen s.t. EM-model predicts sediment load = Omegai
  qsCc = + Ki*eps_s/ws*mean(utot.^3*udelta(1))/(g*(s-1)*(1-psed));    % current-driven above_WBL suspended load

  % slope-driven transport
  Kbdenom1 = + eps_b/tan_phi*mean(utot.^2);
  Kbdenom2 = - slopeFact*eps_b/tan_phi^2*mean(utot.^2)*tanbeta;
  Omega = Omegac*Tc/T + Omegat*Tt/T;  % total sediment load
  K = Omega*(s-1)*g*d50/(Kbdenom1+Kbdenom2+Ksdenom1+Ksdenom2);  % chosen s.t. EM-model predicts total sed load = Omega
  qsCf = - slopeFact*K*(eps_s/ws)^2*mean(utot.^5)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven suspended load
  qbCf = - slopeFact*K*eps_b/tan_phi^2*mean(utot.^3)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven bed load

else  % OPTION-3, above-WBL transport disabled
  eps_s=0.015;
  uwmo=uw/1.4;
  qsVdA = ( qsc + qst )./T.*sqrt((s-1)*g*d50^3)/(1-psed);
  qsCc=0;
  qsCf=0;
  qbCf=0;
end
qs = qsVdA + qsCc + qsCf + qbCf;

% save all relevant variables for passing to TL model.  Note this runs
% faster if the variables are hard-coded, don't use eval
if(nargout>1)
  workspc=struct;
  workspc.Dstar        =Dstar          ;
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
  workspc.h            =h              ;
  workspc.ksd          =ksd            ;
  workspc.ksw          =ksw            ;
  workspc.lambda       =lambda         ;
  workspc.meta         =meta           ;
  workspc.mlambda      =mlambda        ;
  workspc.mu           =mu             ;
  workspc.neta         =neta           ;
  workspc.nlambda      =nlambda        ;
  workspc.param        =param          ;
  workspc.psed         =psed           ;
  workspc.psihat       =psihat         ;
  workspc.psihatc      =psihatc        ;
  workspc.psihatt      =psihatt        ;
  workspc.qs           =qs             ;
  workspc.qsc          =qsc            ;
  workspc.qst          =qst            ;
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
  workspc.wsc          =wsc            ;
  workspc.wst          =wst            ;
  workspc.c1           =c1             ;
  workspc.f25          =f25            ;
  workspc.theta25      =theta25        ;
  workspc.r            =r              ;
  workspc.fws          =fws            ;

  % accounting of suspended seds above WBL
  workspc.eps_s        =eps_s          ;
%   workspc.eps_b        =eps_b          ;
  workspc.uwmo         =uwmo           ;
  workspc.tanbeta      =tanbeta        ;
  if(~isfield(param,'nosusp') || param.nosusp==0)
    workspc.eps_s        =eps_s          ;
    workspc.eps_b        =eps_b          ;
    workspc.uwmo         =uwmo           ;
    workspc.tanbeta      =tanbeta        ;
  end
  if(isfield(param,'Cc'))
    workspc.qs2          =qs2            ;
    workspc.qs3          =qs3            ;
  elseif(~isfield(param,'nosusp') || param.nosusp==0)
    workspc.tan_phi      =tan_phi        ;
    workspc.slopeFact    =slopeFact      ;
    % workspc.bedSuspRatio =bedSuspRatio   ;
    workspc.Lt     =Lt      ;
    workspc.Lc     =Lc      ;
    workspc.wfracc =wfracc  ;
    workspc.wfract =wfract  ;
    workspc.utot   =utot    ;
    workspc.Omegai =Omegai  ;
    workspc.Omega  =Omega   ;
    workspc.Ksdenom1=Ksdenom1 ;
    workspc.Kbdenom1=Kbdenom1 ;
    workspc.Ksdenom2=Ksdenom2 ;
    workspc.Kbdenom2=Kbdenom2 ;
    workspc.Ki      =Ki       ;
    workspc.K       =K        ;
  end
  workspc.qsVdA        =qsVdA          ;
  workspc.qsCc         =qsCc           ;
  workspc.qsCf         =qsCf           ;
  workspc.qbCf         =qbCf           ;

end

% TEST-CODE: override output variable
% eval(['qs = ' outvar ';']);

end  % end of main function
