function [qs,workspc,Omegatc]=qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc)%,outvar)
%
% [qs,workspc]=qtrans_vanderA(x,d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,ws,Aw,Sw,Uw,param)
%
% Calculates transport following van der A (2013)
%
% INPUTS:
%
% d50    : in meters
% d90    : in meters 
% h      : water depth, m
% tanbeta: bottom slope, dimensionless (see calcTanbeta.m)
% Hrms   : rms wave height, m
% kabs   : wavenumber, scalar, rad/m
% omega  : wave frequency, rad/m
% udelta : near-bed mean velocity, m/s
%          Note: wave angle is assumed to be zero in this code.  To
%          incorporate a nonzero wave angle, you should rotate udelta to
%          align with the wave direction.
% ws     : sediment settling velocity, m/s
% Aw,Sw,Uw : wave shape params (as calcualted by Uwave_ruessink2012_params())
% param.{alpha,xi,m,n}: tuning parameters struct
%        alpha: phase lag effect, eqns (24)-(28).  Default 8.2
%        xi   : phase lag effect, eqns (24)-(28).  Default 1.7.  O(1) according to Kranenburg (2013)
%        m    : MPM leading coefficient, eqn (36).  Default 11.0
%        n    : MPM exponent, eqn (36).  Default 1.2
%
%        OPTIONAL: If Cc,Cf are included, above-WBL transport will be
%        calculated explicitly using energetics model.  If not included
%        (recommended), this contribution will be calculated based on
%        boundary layer height plus energetics.
%
%        Cc   : suspended sediment stirring+undertow effect.  Default 0.01
%        Cf   : suspended sediment stirring+slope effect. Default 0.01, but 0.03 may be good
%
%        OPTIONAL: Set param.nosusp=1 to disable above-WBL transport
%        altogether.  This is not recommended as it will cause unreasonable
%        onshore transport in all but the most calm wave conditions.  The
%        option is included for testing purposes only.
%
% param.streamingType : select from 'v' (van der A, 2013) or 'n' (Nielsen,
% 2006).  The Nielsen formulation uses a larger roughness for calculating
% bed streaming, and therefore produces a larger streaming.
%
% OUTPUTS:
%
% qs: total sediment transport flux in m2/s.  A factor 1/(1-psed) is
% included so qs is in terms of bed volume, not sediment volume.
%
% workspc: this is intended only for passing internally to
% {tl,ad}_qtrans_vanderA.m, not for end-user.  Note it is a vector of
% structs, each of which contains bkgd NL variables.
%
%     If you wish to get user access to output variables other than q, it
%     would probably be better to add them to the list of outputs rather
%     than attempt to extract them from the 'workspc' struct.
%

if(isfield(param,'Cc'))
  warning(['Explicit EEM-based suspended-load transport is being depreciated at time of writing.  ' ...
           'Instead, if you remove Cc from params struct, the suspended-load contribution ' ...
           'will be calculated automatically']);
end

% this wrapper loop serves to handle vector inputs
nx=length(h);
for i=1:nx

  % for masked points, make a dummy output with all fields set to 0
  if(Hrms(i)==0)
    blankwksp=workspc(i-1);  % copy as a template.  Can safely assume i>1
    fld=fields(blankwksp);
    for i1=1:length(blankwksp)
      this=getfield(blankwksp,fld{i1});
      if(isstruct(this))
        this=getfield(workspc(i-1),fld{i1});
        fld2=fields(this);
        for i2=1:length(fld2)
          this=setfield(this,fld2{i2},0);
        end
      else
        this=0;
      end
      blankwksp=setfield(blankwksp,fld{i1},this);
    end
    workspc(i)=blankwksp;
    qs(i)=0;
  else
    [qs(i),workspc(i)] = qtrans_vanderA_main(d50(i),d90(i),h(i),tanbeta(i),Hrms(i),kabs(i),omega,udelta(i,:),delta(i),ws(i),Aw(i),Sw(i),Uw(i),param);%,outvar);
  end

end
qs=qs(:);
workspc=workspc(:);

end  % end of wrapper function, start of main function

function [qs,workspc,Omegatc]=qtrans_vanderA_main(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc)%,outvar)

physicalConstants;

% fixed constants
nt=1000;  % for calculation of intra-wave velocity
% delta=0.2;  % OLD: fixed delta=0.2m, per VDA paper

% derived params
Hmo=Hrms*1.4;
c=omega/kabs;

% note, van der A specifies to use significant orbital velocity amplitude,
% while Ruessink et al. (2012) uses rms.  Uwave_ruessink2012() follows the
% Ruessink et al. (2012) convention, so I need to revert to van der A's
% convention here as a special case.
Uw=1.4*Uw;

% intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% relevant to van Der A
t=linspace(0,2*pi/omega,nt);
[uw,uwave_wksp]=Uwave_ruessink2012(omega*t,Aw,Sw,Uw);
r_r2012  =uwave_wksp.r  ;
phi_r2012=uwave_wksp.phi;
uhat=sqrt(2*mean(uw.^2));   % rms wave velocity for full wave cycle, eqn 8
uhatc=max(uw);  % max velocity under crest, see text below eqn 7
uhatt=-min(uw);  % max velocity under trough (positive valued)

% timing of wave velocity direction change, crest, and trough, based on
% Ruessink et al 2012.
phiuc=asin(-r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)))/omega;   % phase at first upcrossing
phidc=pi-phiuc;  % phase at first downcrossing
tuc=phiuc/omega;  % time of first upcrossing
tdc=phidc/omega;  % time of first downcrossing
[~,icu_guess]=max(uw);  % 1st guess
[~,itu_guess]=min(uw);  % 1st guess
tcr=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(icu_guess));  % time of crest
ttr=Uwave_ruessink2012_tcrest(omega,r_r2012,phi_r2012,t(itu_guess));  % time of trough
T=2*pi/omega;  % full wave period
Tc = tdc-tuc;  % duration of crest
Tt=T-Tc;     % duration of trough
Ttu=ttr-tdc;  % duration of deceleration under trough
Tcu=tcr-tuc;  % duration of acceleration under crest

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

% TEST: parameterization of Kranenburg et al. (2012) for boundary layer
% streaming velocity.  Try using this to augment or replace the mean flow at
% top of boundary layer.
% ks=.0082;  % Reniers et al. 2004, calibrated for Duck (Table 4)
% U0 = uhat.^2./c.*( .345 + .7*(ahat/ks).^-.9 - .25./sinh(kabs.*h).^2*0 );
% U0 = 3/4*uhat.^2./c;  % longuet higgins 1958
% udelta(1) = U0 + udelta(1);
% udelta(1)=udelta(1)/4;  % TEST
% udelta(2)=0;  % TEST

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

% % v0: Use linear theory transfer function to get vertical fluid velocities
% worbc=+max(etawc/c*diff(uw)./diff(t));
% worbt=-min(etawt/c*diff(uw)./diff(t));

% Use stokes 2nd order theory to get vertical fluid velocities and correct
% settling velocity.  Follows Malarkey & Davies (2012).  Some of this is
% simply ported from COAWST code.  NOTE, it is a bit odd to do hydrodynamic
% calculations in qtrans_vanderA (by convention, hydro stuff is ideally done
% by callers of this function), but the code can't be easily moved outside
% because it uses ripple height as input.
b=1/r_r2012*(1-sqrt(1-r_r2012^2));  % see MD2012, line after eqn 13b
RR=0.5*(1+b*sin(-phi_r2012));  % see MD2012, eqn 17, note their phi convention is negated
worb1c=pi*Hmo*etawc/(T*h);
worb1t=pi*Hmo*etawt/(T*h);
worb2c=2*worb1c*(2*RR-1);
worb2t=2*worb1t*(2*RR-1);
t1ca = worb1c^2 + 32*worb2c^2;
t1ta = worb1t^2 + 32*worb2t^2;
t1cb = worb1c - sqrt(t1ca);
t1tb = worb1t - sqrt(t1ta);
t1c = 64 + t1cb^2/worb2c^2;
t1t = 64 + t1tb^2/worb2t^2;
t2cb = -worb1c + sqrt( worb1c^2 + 32*worb2c^2 );
t2tb = -worb1t + sqrt( worb1t^2 + 32*worb2t^2 );
t2ca = .125*t2cb/worb2c;
t2ta = .125*t2tb/worb2t;
t2c = 2*acos(t2ca);
t2t = 2*acos(t2ta);
worbc = .125*worb1c*sqrt(t1c) + worb2c*sin(t2c);
worbt = .125*worb1t*sqrt(t1t) + worb2t*sin(t2t);

% apply orbital vels to adjust sediment velocity
wsc=max(ws-worbc,0.001);
wst=max(ws+worbt,0.001);

% sediment load, eqns 23-28.  Note there seems to be a typo in eqns 27-28,
% the bracketed quotient has an issue with dimensions... I changed it into
% what I think is intended.
Pc = param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc);  % eqn 27
Pt = param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst);  % eqn 28
if(abs(thetac)>theta_cr)
  Omegac=param.m*(abs(thetac)-theta_cr).^param.n;  % eqn 2
  %indt=find(uw>0);  % TEST CODE
  %theta_timedep=.5*fwdc.*(uw(indt)+udelta(1)).^2/((s-1)*g*d50);  % TEST CODE
  %Omegac=1/Tc*trapz(t(indt),11*max(0,abs(theta_timedep)-theta_cr).^1.65);  % TEST CODE
else
  Omegac=0;
end
if(abs(thetat)>theta_cr)
  Omegat=param.m*(abs(thetat)-theta_cr).^param.n;  % eqn 2
  %indt=find(uw<0);  % TEST CODE
  %theta_timedep=.5*fwdc.*(uw(indt)+udelta(1)).^2/((s-1)*g*d50);  % TEST CODE
  %Omegat=1/Tt*trapz(t(indt),11*max(0,abs(theta_timedep)-theta_cr).^1.65);  % TEST CODE
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
  qs2=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^3*udelta(1));
  qs3=mean(sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).^5);
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
  if(Omegat>0)
    Lc=max(eps,fzero(@(L)Omegac*d50*exp(-deltasc/L)-0.08*L,.1));
  else
    Lc=0;
  end
  if(Lc>0)
    wfracc = (1-exp(-delta/Lc));
  else
    wfracc=0;
  end
  if(Lt>0)
    wfract = (1-exp(-delta/Lt));
  else
    wfract=0;
  end
  qsVdA = ( qsc*wfracc + qst*wfract )./T.*sqrt((s-1)*g*d50^3)/(1-psed);

  % % OLD: current-driven transport, simply multiply above-WBL sediment load
  % % times udelta.  This kinda works but bars tend to become too peaky.
  % % The slope-driven suspended load transport is missing which would tend to
  % % smooth out the morphology.
  % qsCc = ( (1-wfracc)*Tc./T*Omegac*d50*udelta(1) ...
  %          + (1-wfract)*Tt/T*Omegat*d50*udelta(1) )/(1-psed);  % current-driven
  % qsCf=0;  % slope-driven, left out of this code
  % qs = qs + qsCc + qsCf;  % add suspended contribution to total

  % current- and slope-driven transport are partitioned according to
  % energetics model (i.e., "wave stirring" doing work against gravitational
  % settling).  The two effects are assumed share a common drag coefficient,
  % (K = Cc (==Cf) in Dubarbier's notation).
  eps_s=0.015;
  uwmo=uw/1.4;
  Omegas = (1-wfracc)*Omegac*Tc/T + (1-wfract)*Omegat*Tt/T;  % total load above-WBL
  utot=sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2);  % total current magnitude (time dependent)
  Kdenom1 = + eps_s/ws*mean(utot.^3);
  Kdenom2 = - (eps_s/ws)^2*mean(utot.^4)*tanbeta;
  Kdenom=Kdenom1+Kdenom2;
  K = Omegas*(s-1)*g*d50/(Kdenom1+Kdenom2);  % chosen s.t. EM-model predicts sediment load = Omegas
  qsCc = + K*eps_s/ws*mean(utot.^3*udelta(1))/(g*(s-1)*(1-psed));    % current-driven
  qsCf = - K*(eps_s/ws)^2*mean(utot.^5)*tanbeta/(g*(s-1)*(1-psed));  % slope-driven

else  % OPTION-3, above-WBL transport disabled
  qsVdA = ( qsc + qst )./T.*sqrt((s-1)*g*d50^3)/(1-psed);
  qsCc=0;
  qsCf=0;
end
qs = qsVdA + qsCc + qsCf;

% save all relevant variables for passing to TL model.  Note this runs
% faster if the variables are hard-coded, don't use eval
if(nargout>1)
  workspc=struct;
  workspc.Dstar        =Dstar          ;
  workspc.Hmo          =Hmo            ;
  workspc.Hrms         =Hrms           ;
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
  workspc.kabs         =kabs           ;
  workspc.ksd          =ksd            ;
  workspc.ksw          =ksw            ;
  workspc.lambda       =lambda         ;
  workspc.meta         =meta           ;
  workspc.mlambda      =mlambda        ;
  workspc.mu           =mu             ;
  workspc.neta         =neta           ;
  workspc.nlambda      =nlambda        ;
  workspc.nt           =nt             ;
  workspc.omega        =omega          ;
  workspc.param        =param          ;
  workspc.phi_r2012    =phi_r2012      ;
  workspc.r_r2012      =r_r2012        ;
  workspc.psed         =psed           ;
  workspc.psihat       =psihat         ;
  workspc.psihatc      =psihatc        ;
  workspc.psihatt      =psihatt        ;
  workspc.qs           =qs             ;
  workspc.qsc          =qsc            ;
  workspc.qst          =qst            ;
  workspc.t            =t              ;
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
  workspc.worbc        =worbc          ;
  workspc.worbt        =worbt          ;
  workspc.uwave_wksp   =uwave_wksp     ;
  workspc.Aw           =Aw             ;
  workspc.Sw           =Sw             ;
  workspc.Uw           =Uw             ;
  workspc.c1           =c1             ;
  workspc.b            =b              ;
  workspc.RR           =RR             ;
  workspc.worb1c       =worb1c         ;
  workspc.worb1t       =worb1t         ;
  workspc.worb2c       =worb2c         ;
  workspc.worb2t       =worb2t         ;
  workspc.t1ca         =t1ca           ;
  workspc.t1ta         =t1ta           ;
  workspc.t1cb         =t1cb           ;
  workspc.t1tb         =t1tb           ;
  workspc.t1c          =t1c            ;
  workspc.t1t          =t1t            ;
  workspc.t2cb         =t2cb           ;
  workspc.t2tb         =t2tb           ;
  workspc.t2ca         =t2ca           ;
  workspc.t2ta         =t2ta           ;
  workspc.t2c          =t2c            ;
  workspc.t2t          =t2t            ;
  workspc.worbc        =worbc          ;
  workspc.worbt        =worbt          ;
  workspc.f25          =f25            ;
  workspc.theta25      =theta25        ;
  workspc.r            =r              ;
  workspc.fws          =fws            ;
  workspc.phiuc        =phiuc          ;
  workspc.phidc        =phidc          ;
  workspc.icu_guess    =icu_guess      ;
  workspc.itu_guess    =itu_guess      ;

  % accounting of suspended seds above WBL
  workspc.eps_s        =eps_s          ;
  workspc.uwmo         =uwmo           ;
  workspc.tanbeta      =tanbeta        ;
  if(isfield(param,'Cc'))
    workspc.qs2          =qs2            ;
    workspc.qs3          =qs3            ;
  elseif(~isfield(param,'nosusp') || param.nosusp==0)
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

% TEST-CODE: override output variable
% eval(['qs = ' outvar ';']);

end  % end of main function
