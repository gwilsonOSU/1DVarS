function [qs,workspc]=qtrans_vanderA(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param)
%
% [qs,workspc]=qtrans_vanderA(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param)
%
% Calculates transport following van der A (2013)
%
% The following effects are neglected:
% - vertical velocities (w_min = w_max = 0 in eqns 29-30).  Hence, this is for tunnel flow ONLY.
% - wave angle is not considered in calculating wave-related velocities.
%   See TODO notes in comments, and there may be other places
%
% INPUTS:
%
% d50   : in meters
% d90   : in meters 
% wsc,wst : particle settling velocities under crest and trough, m/s,
%           potentially including a correction for orbital velocity
% udelta: near-bed mean velocity, 2x1 vector, m/s
% uhat=sqrt(2*mean(uw.^2)) : rms wave velocity for full wave cycle, eqn 8, m/s
% uhatc=max(uw) : max velocity under crest, see text below eqn 7, m/s
% uhatt=min(uw) : max velocity under trough (negative valued), m/s
% T   : full wave period, s
% Tc  : duration of crest, s
% Tt  : duration of trough, s
% Ttu : duration of deceleration under trough, s
% Tcu : duration of acceleration under crest, s
% c   : wave phase speed, m/s
% param.{alpha,xi,m,n}: tuning parameters struct
%        alpha: phase lag effect, eqns (24)-(28).  Default 8.2
%        xi   : phase lag effect, eqns (24)-(28).  Default 1.7.  O(1) according to Kranenburg (2013)
%        m    : MPM leading coefficient, eqn (36).  Default 11.0
%        n    : MPM exponent, eqn (36).  Default 1.2
%
% OUTPUTS:
%
% qs: total sediment transport flux in m2/s
%
% workspc: this is intended only for passing internally to
% {tl,ad}_qtrans_vanderA.m, not for end-user.  Note it is a vector of
% structs, each of which contains bkgd NL variables.
%
%     If you wish to get user access to output variables other than q, it
%     would probably be better to add them to the list of outputs rather
%     than attempt to extract them from the 'workspc' struct.
%

% this wrapper loop serves to handle vector inputs.  They can be on a
% spatial grid (hence I use 'nx' below), or they can be a list of individual
% wave strokes.
nx=length(T);
for i=1:nx
  [qs(i),workspc(i)] = qtrans_vanderA_main(d50(i),d90(i),wsc(i),wst(i),udelta(i,:),uhat(i),uhatc(i),uhatt(i),T(i),Tc(i),Tt(i),Ttu(i),Tcu(i),c(i),param);
end
qs=qs(:);

end  % end of wrapper function, start of main function

function [qs,workspc]=qtrans_vanderA_main(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param)

physicalConstants;

% fixed constants
delta=0.2;  % vanderA's approximation, see Fig. 14 and discussion thereof
nt=1000;  % for calculation of intra-wave velocity

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
psihatc=(uhatc)^2/((s-1)*g*d50);  % see Appendix B.  factor 1.27 is for irregular flow
psihatt=(uhatt)^2/((s-1)*g*d50);
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
ucrvec=utildecr*[+1 0]+udelta;  % eqn 12, TODO: wave angle?
utrvec=utildetr*[-1 0]+udelta;  % eqn 13
ucrabs=sqrt(ucrvec(1)^2+ucrvec(2)^2);
utrabs=sqrt(utrvec(1)^2+utrvec(2)^2);
thetac=.5*fwdc.*ucrabs.^2/((s-1)*g*d50);  % eqn 17
thetat=.5*fwdt.*utrabs.^2/((s-1)*g*d50);  % eqn 17
fwd = alpha*fd+(1-alpha)*fw;
alphaw = 4/(3*pi);
thetacx = abs(thetac)*ucrvec(1)/ucrabs;  % eqn 15.  Note, streaming effect excluded for this code
thetatx = abs(thetat)*utrvec(1)/utrabs;  % eqn 15

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

% sediment load, eqns 23-28.  Note there seems to be a typo in eqns 27-28,
% the bracketed quotient has an issue with dimensions... I changed it into
% what I think is intended.
if(eta==0)
  etawc=deltasc;
  etawt=deltast;
else
  etawc=eta;
  etawt=eta;
end
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
  Omegatc=0;
else
  Omegatc=(1-1./Pt)*Omegat;  % eqn 25
end

% transport, eqn 1
qsc=sqrt(abs(thetac)).*Tc.*(Omegacc+Tc./(2*Tcu).*Omegatc).*thetacx./abs(thetac);
qst=sqrt(abs(thetat)).*Tt.*(Omegatt+Tt./(2*Ttu).*Omegact).*thetatx./abs(thetat);
qs = ( qsc + qst )./T.*sqrt((s-1)*g*d50^3);
% qs = qs/(1-psed);  % account for bed porosity

% save all relevant variables for passing to TL model
% TODO: revise this list
if(nargout>1)
  vname={};
  vname{end+1}='Dstar';
  vname{end+1}='Omegac';
  vname{end+1}='Omegacc';
  vname{end+1}='Omegact';
  vname{end+1}='Omegat';
  vname{end+1}='Omegatc';
  vname{end+1}='Omegatt';
  vname{end+1}='Pc';
  vname{end+1}='Pt';
  vname{end+1}='T';
  vname{end+1}='Tc';
  vname{end+1}='Tcu';
  vname{end+1}='Tt';
  vname{end+1}='Ttu';
  vname{end+1}='ahat';
  vname{end+1}='alpha';
  vname{end+1}='branch_A1';
  vname{end+1}='branch_A4';
  vname{end+1}='branch_A5';
  vname{end+1}='d50';
  vname{end+1}='d90';
  vname{end+1}='delta';
  vname{end+1}='deltasc';
  vname{end+1}='deltast';
  vname{end+1}='eta';
  vname{end+1}='etawc';
  vname{end+1}='etawt';
  vname{end+1}='fd';
  vname{end+1}='fw';
  vname{end+1}='fwc';
  vname{end+1}='fwdc';
  vname{end+1}='fwdt';
  vname{end+1}='fwt';
  vname{end+1}='ksd';
  vname{end+1}='ksw';
  vname{end+1}='lambda';
  vname{end+1}='meta';
  vname{end+1}='mlambda';
  vname{end+1}='mu';
  vname{end+1}='neta';
  vname{end+1}='nlambda';
  vname{end+1}='param';
  vname{end+1}='psed';
  vname{end+1}='psihat';
  vname{end+1}='psihatc';
  vname{end+1}='psihatt';
  vname{end+1}='qs';
  vname{end+1}='qsc';
  vname{end+1}='qst';
  vname{end+1}='theta_av';
  vname{end+1}='theta_cr';
  vname{end+1}='thetac';
  vname{end+1}='thetacx';
  vname{end+1}='thetatx';
  vname{end+1}='thetahatc';
  vname{end+1}='thetahatt';
  vname{end+1}='thetat';
  vname{end+1}='ucrabs';
  vname{end+1}='ucrvec';
  vname{end+1}='udabs';
  vname{end+1}='udelta';
  vname{end+1}='uhat';
  vname{end+1}='uhatc';
  vname{end+1}='uhatt';
  vname{end+1}='utildecr';
  vname{end+1}='utildetr';
  vname{end+1}='utrabs';
  vname{end+1}='utrvec';
  vname{end+1}='wsc';
  vname{end+1}='wst';
  vname{end+1}='c';
  vname{end+1}='c1';
  workspc=struct;
  for i=1:length(vname)
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end

end  % end of main function