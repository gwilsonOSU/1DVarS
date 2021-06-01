function [ad_d50,ad_d90,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_param] = ...
    ad_qtrans_vanderA(ad_qs,bkgd,eparam)%,invar)
%
% AD code for qtrans_vanderA.m
%
% OPTIONAL: Set eparam==1 to compute parameter adjoint sensitivity at each
% gridpoint individually.  Else the parameters will be treated as scalar
% constants (i.e., summed to get a single value).
%

if(~exist('eparam'))
  eparam=0;
end

nx=length(ad_qs);

% init ad vars
ad_d50        =zeros(nx,1);
ad_d90        =zeros(nx,1);
ad_ws         =zeros(nx,1);
ad_h          =zeros(nx,1);
ad_Hrms       =zeros(nx,1);
ad_kabs       =zeros(nx,1);
ad_udelta     =zeros(nx,2);
ad_omega=0;
ad_tauwRe=0;
ad_streamingEffect=0;
ad_param.n    =0;
ad_param.m    =0;
ad_param.xi   =0;
ad_param.alpha=0;
ad_param=repmat(ad_param,[nx 1]);

% this wrapper loop serves to handle vector inputs
for i=nx:-1:1
  % tl_qs(i) = ...
  %     tl_qtrans_vanderA(tl_d50,tl_d90,tl_h(i),tl_Hrms(i),tl_kabs(i),...
  %                       tl_udelta(i,:),tl_ws,tl_param,bkgd_qtrans(i));
  [ad1_d50,ad1_d90,ad1_h,ad1_Hrms,ad1_kabs,...
   ad1_omega,ad1_udelta,ad1_ws,ad1_param] = ...
      ad_qtrans_vanderA_main(ad_qs(i),bkgd(i));%,invar);
  ad_d50(i)=ad_d50(i)+ad1_d50   ;
  ad_d90(i)=ad_d90(i)+ad1_d90   ;
  ad_ws(i) =ad_ws(i) +ad1_ws    ;
  ad_h(i)=ad_h(i)+ad1_h;
  ad_Hrms(i)  =ad_Hrms(i)  +ad1_Hrms  ;
  ad_kabs(i)  =ad_kabs(i)  +ad1_kabs  ;
  ad_omega=ad_omega+ad1_omega;
  ad_udelta(i,:)=ad_udelta(i,:)+ad1_udelta;
  ad_param(i).n    =ad_param(i).n    +ad1_param.n    ;
  ad_param(i).m    =ad_param(i).m    +ad1_param.m    ;
  ad_param(i).xi   =ad_param(i).xi   +ad1_param.xi   ;
  ad_param(i).alpha=ad_param(i).alpha+ad1_param.alpha;
end

if(~eparam)
  adp2=struct;
  adp2.n    =sum([ad_param.n    ]);
  adp2.m    =sum([ad_param.m    ]);
  adp2.xi   =sum([ad_param.xi   ]);
  adp2.alpha=sum([ad_param.alpha]);
  ad_param=adp2;
end

end  % end of wrapper function, start of main function

function [ad_d50,ad_d90,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_param]=ad_qtrans_vanderA_main(ad_qs,bkgd)%,invar)

physicalConstants;

% break out NL background vars
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
alpha=bkgd.alpha;
mu=bkgd.mu;

%------------------------------------
% begin AD code
%------------------------------------

% init ad vars
ad_d50   =0;
ad_d90   =0;
ad_h     =0;
ad_Hmo=0;
ad_Hrms  =0;
ad_kabs  =0;
ad_udelta=[0 0];
ad_ws    =0;
ad_param.m    =0;
ad_param.n    =0;
ad_param.alpha=0;
ad_param.xi   =0;
ad_term1=0;
ad_term2=0;
ad_term3=0;
ad_absthetat=0;
ad_absthetac=0;
ad_theta_cr=0;
ad_thetat=0;
ad_thetac=0;
ad_theta_av=0;
ad_thetahatt=0;
ad_thetahatc=0;
ad_Pt=0;
ad_Pc=0;
ad_utrabs=0;
ad_ucrabs=0;
ad_utrvec=[0 0];
ad_ucrvec=[0 0];
ad_utildetr=0;
ad_utildecr=0;
ad_alpha=0;
ad_mu=0;
ad_eta=0;
ad_lambda=0;
ad_mlambda=0;
ad_nlambda=0;
ad_psihat=0;
ad_psihatt=0;
ad_psihatc=0;
ad_Dstar=0;
ad_uw=zeros(1,nt);
ad_meta=0;
ad_neta=0;
ad_fd=0;
ad_fw=0;
ad_ksd=0;
ad_ksw=0;
ad_fwt=0;
ad_fwc=0;
ad_ahat=0;
ad_udabs=0;
ad_uhat=0;
ad_uhatt=0;
ad_uhatc=0;
ad_c=0;
ad_deltast=0;
ad_deltasc=0;
ad_fwdt=0;
ad_fwdc=0;
ad_T=0;
ad_Tt=0;
ad_Tc=0;
ad_Ttu=0;
ad_Tcu=0;
ad_Omegac=0;
ad_Omegat=0;
ad_Omegacc=0;
ad_Omegatt=0;
ad_Omegact=0;
ad_Omegatc=0;
ad_asinarg=0;
ad_uw2mean=0;
ad_r_r2012=0;
ad_phi_r2012=0;
ad_omega=0;
ad_qsc=0;
ad_qst=0;
ad_thetacx=0;
ad_thetatx=0;
ad_etawc=0;
ad_etawt=0;
ad_wsc=0;
ad_wst=0;
ad_worbc=0;
ad_worbt=0;
ad_worb1c=0;
ad_worb1t=0;
ad_worb2c=0;
ad_worb2t=0;
ad_t1c=0;
ad_t1t=0;
ad_t2c=0;
ad_t2t=0;
ad_t1ca=0;
ad_t1ta=0;
ad_t2ca=0;
ad_t2ta=0;
ad_RR=0;
ad_b=0;
ad_streamingEffect=0;
ad_tauwRe=0;
ad_fwd=0;
ad_argc=0;
ad_argt=0;
ad_argc1=0;
ad_argt1=0;
ad_argc2=0;
ad_argt2=0;

% % TEST-CODE: override input variable
% if(~strcmp(invar,'qs'))
%   eval(['ad_' invar '=ad_qs;'])
%   ad_qs=0;
% end

%b14 transport, eqn 1
absthetac=abs(thetac);
absthetat=abs(thetat);
term3 = sqrt((s-1)*g*d50^3)/(1-psed);
%4 tl_qs = (tl_qsc + tl_qst)/T*term3 ...
%         - (qsc + qst)/T^2*term3*tl_T ...
%         + (qsc + qst)/T*tl_term3;
ad_qsc  =ad_qsc  + 1/T*term3            *ad_qs;
ad_qst  =ad_qst  + 1/T*term3            *ad_qs;
ad_T    =ad_T    - (qsc + qst)/T^2*term3*ad_qs;
ad_term3=ad_term3+ (qsc + qst)/T        *ad_qs;
ad_qs=0;
%3 tl_term3 = .5/sqrt((s-1)*g*d50^3)*3*(s-1)*g*d50^2*tl_d50/(1-psed);
ad_d50 = ad_d50 + .5/sqrt((s-1)*g*d50^3)*3*(s-1)*g*d50^2/(1-psed)*ad_term3;
ad_term3=0;
%2 tl_qst = .5/sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat*tl_absthetat ...
%          + sqrt(absthetat)*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat*tl_Tt ...
%          + sqrt(absthetat)*Tt*thetatx/absthetat*( ...
%              + tl_Omegatt ...
%              + tl_Tt/(2*Ttu)*Omegact ...
%              - Tt/(2*Ttu)^2*Omegact*2*tl_Ttu ...
%              + Tt/(2*Ttu)*tl_Omegact ) ...
%          + sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*tl_thetatx/absthetat ...
%          - sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat^2*tl_absthetat;
ad_absthetat=ad_absthetat+ .5/sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat*ad_qst;
ad_Tt       =ad_Tt       + sqrt(absthetat)*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat      *ad_qst;
ad_Omegatt  =ad_Omegatt  + sqrt(absthetat)*Tt*thetatx/absthetat                                *ad_qst;
ad_Tt       =ad_Tt       + sqrt(absthetat)*Tt*thetatx/absthetat/(2*Ttu)*Omegact                *ad_qst;
ad_Ttu      =ad_Ttu      - sqrt(absthetat)*Tt*thetatx/absthetat*Tt/(2*Ttu)^2*Omegact*2         *ad_qst;
ad_Omegact  =ad_Omegact  + sqrt(absthetat)*Tt*thetatx/absthetat*Tt/(2*Ttu)                     *ad_qst;
ad_thetatx  =ad_thetatx  + sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)/absthetat           *ad_qst;
ad_absthetat=ad_absthetat- sqrt(absthetat)*Tt*(Omegatt+Tt/(2*Ttu)*Omegact)*thetatx/absthetat^2 *ad_qst;
ad_qst=0;
%1 tl_qsc = .5/sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac*tl_absthetac ...
%          + sqrt(absthetac)*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac*tl_Tc ...
%          + sqrt(absthetac)*Tc*thetacx/absthetac*( ...
%              + tl_Omegacc ...
%              + tl_Tc/(2*Tcu)*Omegatc ...
%              - Tc/(2*Tcu)^2*Omegatc*2*tl_Tcu ...
%              + Tc/(2*Tcu)*tl_Omegatc ) ...
%          + sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*tl_thetacx/absthetac ...
%          - sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac^2*tl_absthetac;
ad_absthetac=ad_absthetac+ .5/sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac*ad_qsc;
ad_Tc       =ad_Tc       + sqrt(absthetac)*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac      *ad_qsc;
ad_Omegacc  =ad_Omegacc  + sqrt(absthetac)*Tc*thetacx/absthetac                                *ad_qsc;
ad_Tc       =ad_Tc       + sqrt(absthetac)*Tc*thetacx/absthetac/(2*Tcu)*Omegatc                *ad_qsc;
ad_Tcu      =ad_Tcu      - sqrt(absthetac)*Tc*thetacx/absthetac*Tc/(2*Tcu)^2*Omegatc*2         *ad_qsc;
ad_Omegatc  =ad_Omegatc  + sqrt(absthetac)*Tc*thetacx/absthetac*Tc/(2*Tcu)                     *ad_qsc;
ad_thetacx  =ad_thetacx  + sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)/absthetac           *ad_qsc;
ad_absthetac=ad_absthetac- sqrt(absthetac)*Tc*(Omegacc+Tc/(2*Tcu)*Omegatc)*thetacx/absthetac^2 *ad_qsc;
ad_qsc=0;

% %b13 sediment load, eqns 23-28.  Note there seems to be a typo in eqns 27-28,
% % the bracketed quotient has an issue with dimensions... I changed it into
% % what I think is intended
absthetac=abs(thetac);
absthetat=abs(thetat);
if(Pt<=1)
  %   tl_Omegatc=0;
  ad_Omegatc=0;
else
  %   tl_Omegatc = 1./Pt^2*Omegat*tl_Pt ...
  %       + (1-1./Pt)*tl_Omegat;
  ad_Pt    =ad_Pt    + 1./Pt^2*Omegat*ad_Omegatc;
  ad_Omegat=ad_Omegat+ (1-1./Pt)     *ad_Omegatc;
  ad_Omegatc=0;
end
if(Pc<=1)
  %   tl_Omegact=0;
  ad_Omegact=0;
else
  %   tl_Omegact = 1./Pc^2*Omegac*tl_Pc ...
  %       + (1-1./Pc)*tl_Omegac;
  ad_Pc    =ad_Pc    + 1./Pc^2*Omegac*ad_Omegact;
  ad_Omegac=ad_Omegac+ (1-1./Pc)     *ad_Omegact;
  ad_Omegact=0;
end
if(Pt>1)
  %   tl_Omegatt = tl_Omegat./Pt ...
  %       - Omegat./Pt^2*tl_Pt;
  ad_Omegat=ad_Omegat+ 1./Pt       *ad_Omegatt;
  ad_Pt    =ad_Pt    - Omegat./Pt^2*ad_Omegatt;
  ad_Omegatt=0;
else
  %   tl_Omegatt=tl_Omegat;
  ad_Omegat=ad_Omegat+ad_Omegatt;
  ad_Omegatt=0;
end
if(Pc>1)
  %   tl_Omegacc = tl_Omegac./Pc ...
  %       - Omegac./Pc^2*tl_Pc;
  ad_Omegac=ad_Omegac+ 1./Pc       *ad_Omegacc;
  ad_Pc    =ad_Pc    - Omegac./Pc^2*ad_Omegacc;
  ad_Omegacc=0;
else
  %   tl_Omegacc=tl_Omegac;
  ad_Omegac=ad_Omegac+ad_Omegacc;
  ad_Omegacc=0;
end
if(abs(thetat)>theta_cr)
  %   tl_Omegat = tl_param.m*(absthetat-theta_cr).^param.n ...
  %       + param.n*param.m*(absthetat-theta_cr).^(param.n-1)*(tl_absthetat-tl_theta_cr) ...
  %       + param.m*(absthetat-theta_cr).^param.n*log(absthetat-theta_cr)*tl_param.n;
  ad_param.m  =ad_param.m  + (absthetat-theta_cr).^param.n                                *ad_Omegat;
  ad_absthetat=ad_absthetat+ param.n*param.m*(absthetat-theta_cr).^(param.n-1)            *ad_Omegat;
  ad_theta_cr =ad_theta_cr - param.n*param.m*(absthetat-theta_cr).^(param.n-1)            *ad_Omegat;
  ad_param.n  =ad_param.n  + param.m*(absthetat-theta_cr).^param.n*log(absthetat-theta_cr)*ad_Omegat;
  ad_Omegat=0;
else
  %   tl_Omegat=0;
  ad_Omegat=0;
end
if(abs(thetac)>theta_cr)
  %   tl_Omegac = tl_param.m*(absthetac-theta_cr).^param.n ...
  %       + param.n*param.m*(absthetac-theta_cr).^(param.n-1)*(tl_absthetac-tl_theta_cr) ...
  %       + param.m*(absthetac-theta_cr).^param.n*log(absthetac-theta_cr)*tl_param.n;
  ad_param.m  =ad_param.m  + (absthetac-theta_cr).^param.n                                *ad_Omegac;
  ad_absthetac=ad_absthetac+ param.n*param.m*(absthetac-theta_cr).^(param.n-1)            *ad_Omegac;
  ad_theta_cr =ad_theta_cr - param.n*param.m*(absthetac-theta_cr).^(param.n-1)            *ad_Omegac;
  ad_param.n  =ad_param.n  + param.m*(absthetac-theta_cr).^param.n*log(absthetac-theta_cr)*ad_Omegac;
  ad_Omegac=0;
else
  %   tl_Omegac=0;
  ad_Omegac=0;
end
%4 tl_absthetat = sign(thetat)*tl_thetat;
ad_thetat=ad_thetat+sign(thetat)*ad_absthetat;
ad_absthetat=0;
%3 tl_absthetac = sign(thetac)*tl_thetac;
ad_thetac =ad_thetac+ sign(thetac)*ad_absthetac;
ad_absthetac=0;
%2 tl_Pt = tl_param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst) ...
%         + param.alpha*( ...
%             + tl_param.xi*uhatt./c ...
%             + param.xi*tl_uhatt./c ...
%             - param.xi*uhatt./c^2*tl_c ...
%             ).*etawt./(2*(Tt-Ttu)*wst) ...
%         + param.alpha*(1+param.xi*uhatt./c).*tl_etawt./(2*(Tt-Ttu)*wst) ...
%         - param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)^2*( ...
%             + 2*(tl_Tt-tl_Ttu)*wst ...
%             + 2*(Tt-Ttu)*tl_wst );
ad_param.alpha=ad_param.alpha+ (1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)                         *ad_Pt;
ad_param.xi   =ad_param.xi   + param.alpha*etawt./(2*(Tt-Ttu)*wst)*uhatt./c                           *ad_Pt;
ad_uhatt      =ad_uhatt      + param.alpha*etawt./(2*(Tt-Ttu)*wst)*param.xi./c                        *ad_Pt;
ad_c          =ad_c          - param.alpha*etawt./(2*(Tt-Ttu)*wst)*param.xi*uhatt./c^2                *ad_Pt;
ad_etawt      =ad_etawt      + param.alpha*(1+param.xi*uhatt./c)./(2*(Tt-Ttu)*wst)                    *ad_Pt;
ad_Tt         =ad_Tt         - param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)^2*2*wst     *ad_Pt;
ad_Ttu        =ad_Ttu        + param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)^2*2*wst     *ad_Pt;
ad_wst        =ad_wst        - param.alpha*(1+param.xi*uhatt./c).*etawt./(2*(Tt-Ttu)*wst)^2*2*(Tt-Ttu)*ad_Pt;
ad_Pt=0;
%1 tl_Pc = tl_param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc) ...
%         + param.alpha*( ...
%             - tl_param.xi*uhatc./c ...
%             - param.xi*tl_uhatc./c ...
%             + param.xi*uhatc./c^2*tl_c ...
%             ).*etawc./(2*(Tc-Tcu)*wsc) ...
%         + param.alpha*(1-param.xi*uhatc./c).*tl_etawc./(2*(Tc-Tcu)*wsc) ...
%         - param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)^2*( ...
%             + 2*(tl_Tc-tl_Tcu)*wsc ...
%             + 2*(Tc-Tcu)*tl_wsc );
ad_param.alpha=ad_param.alpha+ (1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)                         *ad_Pc;
ad_param.xi   =ad_param.xi   - param.alpha*etawc./(2*(Tc-Tcu)*wsc)*uhatc./c                           *ad_Pc;
ad_uhatc      =ad_uhatc      - param.alpha*etawc./(2*(Tc-Tcu)*wsc)*param.xi./c                        *ad_Pc;
ad_c          =ad_c          + param.alpha*etawc./(2*(Tc-Tcu)*wsc)*param.xi*uhatc./c^2                *ad_Pc;
ad_etawc      =ad_etawc      + param.alpha*(1-param.xi*uhatc./c)./(2*(Tc-Tcu)*wsc)                    *ad_Pc;
ad_Tc         =ad_Tc         - param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)^2*2*wsc     *ad_Pc;
ad_Tcu        =ad_Tcu        + param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)^2*2*wsc     *ad_Pc;
ad_wsc        =ad_wsc        - param.alpha*(1-param.xi*uhatc./c).*etawc./(2*(Tc-Tcu)*wsc)^2*2*(Tc-Tcu)*ad_Pc;
ad_Pc=0;

% %b13a Use stokes 2nd order theory to get vertical fluid velocities and correct
% % settling velocity.  Follows Malarkey & Davies (2012).  Some of this is
% % simply ported from COAWST code.
b=1/r_r2012*(1-sqrt(1-r_r2012^2));  % see MD2012, line after eqn 13b
RR=0.5*(1+b*sin(-phi_r2012));  % see MD2012, eqn 17, note their phi convention is negated
worb1c=pi*Hmo*etawc/(T*h);
worb1t=pi*Hmo*etawt/(T*h);
worb2c=2*worb1c*(2*RR-1);
worb2t=2*worb1t*(2*RR-1);
t1ca=-worb1c+sqrt(worb1c^2+32*worb2c^2);
t1c=64-t1ca^2/worb2c^2;
t2ca=t1ca/worb2c;
t2c=acos(.125*t2ca);
worbc = .125*worb1c*sqrt(t1c) ...
       + worb2c*sin(2*t2c);
t1ta=-worb1t+sqrt(worb1t^2+32*worb2t^2);
t1t=64-t1ta^2/worb2t^2;
t2ta=t1ta/worb2t;
t2t=acos(.125*t2ta);
worbt = .125*worb1t*sqrt(t1t) ...
       + worb2t*sin(2*t2t);
wsc=ws-worbc;
wst=max(ws+worbt,0);

%18 tl_wst = tl_ws + tl_worbt;
ad_ws   =ad_ws   + ad_wst;
ad_worbt=ad_worbt+ ad_wst;
ad_wst=0;
%17 tl_wsc = tl_ws - tl_worbc;
ad_ws   =ad_ws   + ad_wsc;
ad_worbc=ad_worbc- ad_wsc;
ad_wsc=0;
%16 tl_worbt = .125*sqrt(t1t)*tl_worb1t ...
%     + .5*.125*worb1t/sqrt(t1t)*tl_t1t ...
%     + sin(2*t2t)*tl_worb2t ...
%     + worb2t*cos(2*t2t)*2*tl_t2t;
ad_worb1t=ad_worb1t+ .125*sqrt(t1t)          *ad_worbt;
ad_t1t   =ad_t1t   + .5*.125*worb1t/sqrt(t1t)*ad_worbt;
ad_worb2t=ad_worb2t+ sin(2*t2t)              *ad_worbt;
ad_t2t   =ad_t2t   + worb2t*cos(2*t2t)*2     *ad_worbt;
ad_worbt=0;
%15 tl_t2t = -1/sqrt(1-(.125*t2ta)^2)*.125*tl_t2ta;
ad_t2ta=ad_t2ta-1/sqrt(1-(.125*t2ta)^2)*.125*ad_t2t;
ad_t2t=0;
%14 tl_t2ta = tl_t1ta/worb2t ...
%           - t1ta/worb2t^2*tl_worb2t;
ad_t1ta  =ad_t1ta  + 1/worb2t     *ad_t2ta;
ad_worb2t=ad_worb2t- t1ta/worb2t^2*ad_t2ta;
ad_t2ta=0;
%13 tl_t1t = -2*t1ta/worb2t^2*tl_t1ta ...
%          + 2*t1ta^2/worb2t^3*tl_worb2t;
ad_t1ta  =ad_t1ta  - 2*t1ta/worb2t^2  *ad_t1t;
ad_worb2t=ad_worb2t+ 2*t1ta^2/worb2t^3*ad_t1t;
ad_t1t=0;
%12 tl_t1ta = -tl_worb1t ...
%           + .5/sqrt(worb1t^2+32*worb2t^2)*( ...
%               + 2*worb1t*tl_worb1t ...
%               + 2*32*worb2t*tl_worb2t );
ad_worb1t=ad_worb1t+ (-1+.5/sqrt(worb1t^2+32*worb2t^2)*2*worb1t)*ad_t1ta;
ad_worb2t=ad_worb2t+ .5/sqrt(worb1t^2+32*worb2t^2)*2*32*worb2t  *ad_t1ta;
ad_t1ta=0;
%11 tl_worbc = .125*sqrt(t1c)*tl_worb1c ...
%     + .5*.125*worb1c/sqrt(t1c)*tl_t1c ...
%     + sin(2*t2c)*tl_worb2c ...
%     + worb2c*cos(2*t2c)*2*tl_t2c;
ad_worb1c=ad_worb1c+ .125*sqrt(t1c)          *ad_worbc;
ad_t1c   =ad_t1c   + .5*.125*worb1c/sqrt(t1c)*ad_worbc;
ad_worb2c=ad_worb2c+ sin(2*t2c)              *ad_worbc;
ad_t2c   =ad_t2c   + worb2c*cos(2*t2c)*2     *ad_worbc;
ad_worbc=0;
%10 tl_t2c = -1/sqrt(1-(.125*t2ca)^2)*.125*tl_t2ca;
ad_t2ca=ad_t2ca-1/sqrt(1-(.125*t2ca)^2)*.125*ad_t2c;
ad_t2c=0;
%9 tl_t2ca = tl_t1ca/worb2c ...
%           - t1ca/worb2c^2*tl_worb2c;
ad_t1ca  =ad_t1ca  + 1/worb2c     *ad_t2ca;
ad_worb2c=ad_worb2c- t1ca/worb2c^2*ad_t2ca;
ad_t2ca=0;
%8 tl_t1c = -2*t1ca/worb2c^2*tl_t1ca ...
%          + 2*t1ca^2/worb2c^3*tl_worb2c;
ad_t1ca  =ad_t1ca  - 2*t1ca/worb2c^2  *ad_t1c;
ad_worb2c=ad_worb2c+ 2*t1ca^2/worb2c^3*ad_t1c;
ad_t1c=0;
%7 tl_t1ca = -tl_worb1c ...
%           + .5/sqrt(worb1c^2+32*worb2c^2)*( ...
%               + 2*worb1c*tl_worb1c ...
%               + 2*32*worb2c*tl_worb2c );
ad_worb1c=ad_worb1c+ (-1+.5/sqrt(worb1c^2+32*worb2c^2)*2*worb1c)*ad_t1ca;
ad_worb2c=ad_worb2c+ .5/sqrt(worb1c^2+32*worb2c^2)*2*32*worb2c  *ad_t1ca;
ad_t1ca=0;
%6 tl_worb2t = 2*(2*RR-1)*tl_worb1t ...
%     + 2*worb1t*2*tl_RR;
ad_worb1t=ad_worb1t+ 2*(2*RR-1)*ad_worb2t;
ad_RR    =ad_RR    + 2*worb1t*2*ad_worb2t;
ad_worb2t=0;
%5 tl_worb2c = 2*(2*RR-1)*tl_worb1c ...
%     + 2*worb1c*2*tl_RR;
ad_worb1c=ad_worb1c+ 2*(2*RR-1)*ad_worb2c;
ad_RR    =ad_RR    + 2*worb1c*2*ad_worb2c;
ad_worb2c=0;
%4 tl_worb1t = pi*etawt/(T*h)*tl_Hmo ...
%     + pi*Hmo/(T*h)*tl_etawt ...
%     - pi*Hmo*etawt/(T^2*h)*tl_T ...
%     - pi*Hmo*etawt/(T*h^2)*tl_h;
ad_Hmo  =ad_Hmo  + pi*etawt/(T*h)      *ad_worb1t;
ad_etawt=ad_etawt+ pi*Hmo/(T*h)        *ad_worb1t;
ad_T    =ad_T    - pi*Hmo*etawt/(T^2*h)*ad_worb1t;
ad_h    =ad_h    - pi*Hmo*etawt/(T*h^2)*ad_worb1t;
ad_worb1t=0;
%3 tl_worb1c = pi*etawc/(T*h)*tl_Hmo ...
%     + pi*Hmo/(T*h)*tl_etawc ...
%     - pi*Hmo*etawc/(T^2*h)*tl_T ...
%     - pi*Hmo*etawc/(T*h^2)*tl_h;
ad_Hmo  =ad_Hmo  + pi*etawc/(T*h)      *ad_worb1c;
ad_etawc=ad_etawc+ pi*Hmo/(T*h)        *ad_worb1c;
ad_T    =ad_T    - pi*Hmo*etawc/(T^2*h)*ad_worb1c;
ad_h    =ad_h    - pi*Hmo*etawc/(T*h^2)*ad_worb1c;
ad_worb1c=0;
%2 tl_RR = 0.5*sin(-phi_r2012)*tl_b ...
%         - 0.5*b*cos(-phi_r2012)*tl_phi_r2012;
ad_b        =ad_b        + 0.5*sin(-phi_r2012)  *ad_RR;
ad_phi_r2012=ad_phi_r2012- 0.5*b*cos(-phi_r2012)*ad_RR;
ad_RR=0;
%1 tl_b = -1/r_r2012^2*(1-sqrt(1-r_r2012^2))*tl_r_r2012 ...
%        + .5/r_r2012/sqrt(1-r_r2012^2)*2*r_r2012*tl_r_r2012;
ad_r_r2012 = ad_r_r2012 + (-1/r_r2012^2*(1-sqrt(1-r_r2012^2)) + .5/r_r2012/sqrt(1-r_r2012^2)*2*r_r2012)*ad_b;
ad_b=0;
if(eta==0)
  % tl_etawt=tl_deltast;
  ad_deltast=ad_deltast+ad_etawt;
  ad_etawt=0;
  % tl_etawc=tl_deltasc;
  ad_deltasc=ad_deltasc+ad_etawc;
  ad_etawc=0;
else
  % tl_etawt=tl_eta;
  ad_eta=ad_eta+ad_etawt;
  ad_etawt=0;
  % tl_etawc=tl_eta;
  ad_eta=ad_eta+ad_etawc;
  ad_etawc=0;
end

% %b12 sheet flow layer thickness, Appendix C
if(d50<=.15e-3)
  %   tl_deltast=tl_d50*25*thetahatt ...
  %       + d50*25*tl_thetahatt;
  ad_d50      =ad_d50      + 25*thetahatt*ad_deltast;
  ad_thetahatt=ad_thetahatt+ d50*25      *ad_deltast;
  ad_deltast=0;
  %   tl_deltasc=tl_d50*25*thetahatc ...
  %       + d50*25*tl_thetahatc;
  ad_d50      =ad_d50      + 25*thetahatc*ad_deltasc;
  ad_thetahatc=ad_thetahatc+ d50*25      *ad_deltasc;
  ad_deltasc=0;
elseif(.15e-3<d50&d50<.2e-3)
  %   tl_deltast = tl_d50*(25-12*(d50-.15e-3)/.05e-3) ...
  %       - d50*12*tl_d50/.05e-3;
  ad_d50=ad_d50+ (25-12*(d50-.15e-3)/.05e-3)*ad_deltast;
  ad_d50=ad_d50- d50*12/.05e-3              *ad_deltast;
  ad_deltast=0;
  %   tl_deltasc = tl_d50*(25-12*(d50-.15e-3)/.05e-3) ...
  %       - d50*12*tl_d50/.05e-3;
  ad_d50=ad_d50+ (25-12*(d50-.15e-3)/.05e-3)*ad_deltasc;
  ad_d50=ad_d50- d50*12/.05e-3              *ad_deltasc;
  ad_deltasc=0;
else
  %   tl_deltast = tl_d50*13*thetahatt ...
  %       + d50*13*tl_thetahatt;
  ad_d50      =ad_d50      + 13*thetahatt*ad_deltast;
  ad_thetahatt=ad_thetahatt+ d50*13      *ad_deltast;
  ad_deltast=0;
  %   tl_deltasc = tl_d50*13*thetahatc ...
  %       + d50*13*tl_thetahatc;
  ad_d50      =ad_d50      + 13*thetahatc*ad_deltasc;
  ad_thetahatc=ad_thetahatc+ d50*13      *ad_deltasc;
  ad_deltasc=0;
end
%2 tl_thetahatt = .5/((s-1)*g)*( tl_fwdt.*uhatt.^2/d50 ...
%                            + 2*fwdt.*uhatt*tl_uhatt/d50 ...
%                            - fwdt.*uhatt.^2/d50^2*tl_d50 );
ad_fwdt =ad_fwdt + .5/((s-1)*g).*uhatt.^2/d50       *ad_thetahatt;
ad_uhatt=ad_uhatt+ .5/((s-1)*g)*2*fwdt.*uhatt/d50   *ad_thetahatt;
ad_d50  =ad_d50  - .5/((s-1)*g)*fwdt.*uhatt.^2/d50^2*ad_thetahatt;
ad_thetahatt=0;
%1 tl_thetahatc = .5/((s-1)*g)*( tl_fwdc.*uhatc.^2/d50 ...
%                            + 2*fwdc.*uhatc*tl_uhatc/d50 ...
%                            - fwdc.*uhatc.^2/d50^2*tl_d50 );
ad_fwdc =ad_fwdc + .5/((s-1)*g).*uhatc.^2/d50       *ad_thetahatc;
ad_uhatc=ad_uhatc+ .5/((s-1)*g)*2*fwdc.*uhatc/d50   *ad_thetahatc;
ad_d50  =ad_d50  - .5/((s-1)*g)*fwdc.*uhatc.^2/d50^2*ad_thetahatc;
ad_thetahatc=0;

% %b11 continued
%5 tl_thetatx = tl_thetat*utrvec(1)/utrabs ...
%     + thetat*tl_utrvec(1)/utrabs ...
%     - thetat*utrvec(1)/utrabs^2*tl_utrabs ...
%     + tl_streamingEffect;
ad_thetat         =ad_thetat         + utrvec(1)/utrabs         *ad_thetatx;
ad_utrvec(1)      =ad_utrvec(1)      + thetat/utrabs            *ad_thetatx;
ad_utrabs         =ad_utrabs         - thetat*utrvec(1)/utrabs^2*ad_thetatx;
ad_streamingEffect=ad_streamingEffect+                         1*ad_thetatx;
ad_thetatx=0;
%4 tl_thetacx = tl_thetac*ucrvec(1)/ucrabs ...
%     + thetac*tl_ucrvec(1)/ucrabs ...
%     - thetac*ucrvec(1)/ucrabs^2*tl_ucrabs ...
%     + tl_streamingEffect;
ad_thetac         =ad_thetac         + ucrvec(1)/ucrabs         *ad_thetacx;
ad_ucrvec(1)      =ad_ucrvec(1)      + thetac/ucrabs            *ad_thetacx;
ad_ucrabs         =ad_ucrabs         - thetac*ucrvec(1)/ucrabs^2*ad_thetacx;
ad_streamingEffect=ad_streamingEffect+                         1*ad_thetacx;
ad_thetacx=0;
%3 tl_streamingEffect = 1/((s-1)*g)*( tl_tauwRe/d50 - tauwRe/d50^2*tl_d50 );
ad_tauwRe=ad_tauwRe+ 1/((s-1)*g)/d50         *ad_streamingEffect;
ad_d50   =ad_d50   - 1/((s-1)*g)*tauwRe/d50^2*ad_streamingEffect;
ad_streamingEffect=0;
%2 tl_tauwRe = .5*alphaw*( ...
%     tl_fwd*uhat^3/c ...
%     + 3*fwd*uhat^2/c*tl_uhat ...
%     - fwd*uhat^3/c^2*tl_c );
ad_fwd =ad_fwd + .5*alphaw*uhat^3/c      *ad_tauwRe;
ad_uhat=ad_uhat+ .5*alphaw*3*fwd*uhat^2/c*ad_tauwRe;
ad_c   =ad_c   - .5*alphaw*fwd*uhat^3/c^2*ad_tauwRe;
ad_tauwRe=0;
%1 tl_fwd = tl_alpha*fd ...
%          + alpha*tl_fd ...
%          - tl_alpha*fw ...
%          + (1-alpha)*tl_fw;
ad_alpha=ad_alpha+ (fd-fw)  *ad_fwd;
ad_fd   =ad_fd   + alpha    *ad_fwd;
ad_fw   =ad_fw   + (1-alpha)*ad_fwd;
ad_fwd=0;

% %b11 other BBL derived parameters
%10 tl_thetat = .5/((s-1)*g)*( tl_fwdt.*utrabs.^2/d50 ...
%                            + 2*fwdt.*utrabs*tl_utrabs/d50 ...
%                            - fwdt.*utrabs.^2/d50^2*tl_d50 );
ad_fwdt  =ad_fwdt  + .5/((s-1)*g).*utrabs.^2/d50       *ad_thetat;
ad_utrabs=ad_utrabs+ .5/((s-1)*g)*2*fwdt.*utrabs/d50   *ad_thetat;
ad_d50   =ad_d50   - .5/((s-1)*g)*fwdt.*utrabs.^2/d50^2*ad_thetat;
ad_thetat=0;
%9 tl_thetac = .5/((s-1)*g)*( tl_fwdc.*ucrabs.^2/d50 ...
%                            + 2*fwdc.*ucrabs*tl_ucrabs/d50 ...
%                            - fwdc.*ucrabs.^2/d50^2*tl_d50 );
ad_fwdc  =ad_fwdc  + .5/((s-1)*g).*ucrabs.^2/d50       *ad_thetac;
ad_ucrabs=ad_ucrabs+ .5/((s-1)*g)*2*fwdc.*ucrabs/d50   *ad_thetac;
ad_d50   =ad_d50   - .5/((s-1)*g)*fwdc.*ucrabs.^2/d50^2*ad_thetac;
ad_thetac=0;
%8 tl_utrabs = .5/sqrt(utrvec(1)^2+utrvec(2)^2)*( ...
%     2*utrvec(1)*tl_utrvec(1) ...
%     + 2*utrvec(2)*tl_utrvec(2) );
coef=.5/sqrt(utrvec(1)^2+utrvec(2)^2);
ad_utrvec(1)=ad_utrvec(1)+ coef*2*utrvec(1)*ad_utrabs;
ad_utrvec(2)=ad_utrvec(2)+ coef*2*utrvec(2)*ad_utrabs;
ad_utrabs=0;
%7 tl_ucrabs = .5/sqrt(ucrvec(1)^2+ucrvec(2)^2)*( ...
%     2*ucrvec(1)*tl_ucrvec(1) ...
%     + 2*ucrvec(2)*tl_ucrvec(2) );
coef=.5/sqrt(ucrvec(1)^2+ucrvec(2)^2);
ad_ucrvec(1)=ad_ucrvec(1)+ coef*2*ucrvec(1)*ad_ucrabs;
ad_ucrvec(2)=ad_ucrvec(2)+ coef*2*ucrvec(2)*ad_ucrabs;
ad_ucrabs=0;
%5 tl_ucrvec(2) = tl_udelta(2);
ad_udelta(2) = ad_udelta(2) + ad_ucrvec(2);
ad_ucrvec(2)=0;
%5 tl_utrvec(2) = tl_udelta(2);
ad_udelta(2) = ad_udelta(2) + ad_utrvec(2);
ad_utrvec(2)=0;
%5 tl_ucrvec(1) = tl_utildecr(1) + tl_udelta(1);
ad_utildecr(1)=ad_utildecr(1)+ ad_ucrvec(1);
ad_udelta(1)  =ad_udelta(1)  + ad_ucrvec(1);
ad_ucrvec(1)=0;
%5 tl_utrvec(1) = tl_utildetr(1) + tl_udelta(1);
ad_utildetr(1)=ad_utildetr(1)+ ad_utrvec(1);
ad_udelta(1)  =ad_udelta(1)  + ad_utrvec(1);
ad_utrvec(1)=0;
%4 tl_fwdt = tl_alpha*fd ...
%        + alpha*tl_fd ...
%        - tl_alpha*fwt ...
%        + (1-alpha)*tl_fwt;
ad_alpha=ad_alpha+ fd       *ad_fwdt;
ad_fd   =ad_fd   + alpha    *ad_fwdt;
ad_alpha=ad_alpha- fwt      *ad_fwdt;
ad_fwt  =ad_fwt  + (1-alpha)*ad_fwdt;
ad_fwdt=0;
%3 tl_fwdc = tl_alpha*fd ...
%           + alpha*tl_fd ...
%           - tl_alpha*fwc ...
%           + (1-alpha)*tl_fwc;
ad_alpha=ad_alpha+ fd       *ad_fwdc;
ad_fd   =ad_fd   + alpha    *ad_fwdc;
ad_alpha=ad_alpha- fwc      *ad_fwdc;
ad_fwc  =ad_fwc  + (1-alpha)*ad_fwdc;
ad_fwdc=0;
if(ahat/ksw>1.587)  % eqn 21

  argc2=2*Tcu/Tc;
  argt2=2*Ttu/Tt;
  argc1=argc2^c1*ahat/ksw;
  argt1=argt2^c1*ahat/ksw;
  argc=5.21*argc1^(-.19);
  argt=5.21*argt1^(-.19);

  %8 tl_fwt=.00251*exp(argt)*tl_argt;
  ad_argt=ad_argt+.00251*exp(argt)*ad_fwt;
  ad_fwt=0;
  %7 tl_fwc=.00251*exp(argc)*tl_argc;
  ad_argc=ad_argc+.00251*exp(argc)*ad_fwc;
  ad_fwc=0;
  %6 tl_argt = -.19*5.21*argt1^(-1.19)*tl_argt1;
  ad_argt1=ad_argt1-.19*5.21*argt1^(-1.19)*ad_argt;
  ad_argt=0;
  %5 tl_argc = -.19*5.21*argc1^(-1.19)*tl_argc1;
  ad_argc1=ad_argc1-.19*5.21*argc1^(-1.19)*ad_argc;
  ad_argc=0;
  %4 tl_argt1 = c1*argt2^(c1-1)*ahat/ksw*tl_argt2 ...
  %     + argt2^c1*tl_ahat/ksw ...
  %     - argt2^c1*ahat/ksw^2*tl_ksw;
  ad_argt2=ad_argt2+ c1*argt2^(c1-1)*ahat/ksw*ad_argt1;
  ad_ahat =ad_ahat + argt2^c1/ksw            *ad_argt1;
  ad_ksw  =ad_ksw  - argt2^c1*ahat/ksw^2     *ad_argt1;
  ad_argt1=0;
  %3 tl_argc1 = c1*argc2^(c1-1)*ahat/ksw*tl_argc2 ...
  %     + argc2^c1*tl_ahat/ksw ...
  %     - argc2^c1*ahat/ksw^2*tl_ksw;
  ad_argc2=ad_argc2+ c1*argc2^(c1-1)*ahat/ksw*ad_argc1;
  ad_ahat =ad_ahat + argc2^c1/ksw            *ad_argc1;
  ad_ksw  =ad_ksw  - argc2^c1*ahat/ksw^2     *ad_argc1;
  ad_argc1=0;
  %2 tl_argt2 = 2*tl_Ttu/Tt - 2*Ttu/Tt^2*tl_Tt;
  ad_Ttu=ad_Ttu+ 2/Tt      *ad_argt2;
  ad_Tt =ad_Tt - 2*Ttu/Tt^2*ad_argt2;
  ad_argt2=0;
  %1 tl_argc2 = 2*tl_Tcu/Tc - 2*Tcu/Tc^2*tl_Tc;
  ad_Tcu=ad_Tcu+ 2/Tc      *ad_argc2;
  ad_Tc =ad_Tc - 2*Tcu/Tc^2*ad_argc2;
  ad_argc2=0;

else
  %   tl_fwc=0;
  %   tl_fwt=0;
  ad_fwc=0;
  ad_fwt=0;
end
%1 tl_alpha = tl_udabs./(udabs+uhat) ...
%     - udabs./(udabs+uhat)^2*( tl_udabs + tl_uhat );
ad_udabs=ad_udabs+ 1./(udabs+uhat)      *ad_alpha;
ad_udabs=ad_udabs- udabs./(udabs+uhat)^2*ad_alpha;
ad_uhat =ad_uhat - udabs./(udabs+uhat)^2*ad_alpha;
ad_alpha=0;

% %b10 shields parameter related parameters.  Requires solving a 5-eqn nonlinear
% % system, so this is done in its own code.
%2 [tl_theta_av,tl_ksd,tl_ksw,tl_fd,tl_fw,A] = ...
%     tl_vanderA_shields(tl_d50,tl_d90,tl_udabs,tl_uhat,...
%                        tl_mu,tl_eta,tl_lambda,tl_ahat,...
%                        d50,d90,udabs,uhat,delta,mu,eta,lambda,ahat,...
%                        ksd,ksw,fd,fw,...
%                        branch_A1,branch_A4,branch_A5);
[ad1_d50,ad1_d90,ad1_udabs,ad1_uhat,...
 ad1_mu,ad1_eta,ad1_lambda,ad1_ahat] = ...
    ad_vanderA_shields(ad_theta_av,ad_ksd,ad_ksw,ad_fd,ad_fw,...
                       d50,udabs,uhat,delta,mu,eta,lambda,...
                       ahat,ksd,ksw,fd,fw,...
                       branch_A1,branch_A4,branch_A5);
ad_d50   =ad_d50   +ad1_d50  ;
ad_d90   =ad_d90   +ad1_d90  ;
ad_udabs=ad_udabs+ad1_udabs;
ad_uhat  =ad_uhat  +ad1_uhat ;
ad_mu    =ad_mu    +ad1_mu   ;
ad_eta   =ad_eta   +ad1_eta  ;
ad_lambda=ad_lambda+ad1_lambda;
ad_ahat  =ad_ahat  +ad1_ahat ;
ad_theta_av=0;
ad_ksd     =0;
ad_ksw     =0;
ad_fd      =0;
ad_fw      =0;

if(sqrt(udelta(1)^2+udelta(2)^2)==0)
  % tl_udabs=0;
  ad_udabs=0;
else
  %1 tl_udabs=.5/sqrt(udelta(1)^2+udelta(2)^2)*( ...
  %     2*udelta(1)*tl_udelta(1) ...
  %     + 2*udelta(2)*tl_udelta(2) );
  coef=.5/sqrt(udelta(1)^2+udelta(2)^2);
  ad_udelta(1)=ad_udelta(1)+ coef*2*udelta(1)*ad_udabs;
  ad_udelta(2)=ad_udelta(2)+ coef*2*udelta(2)*ad_udabs;
  ad_udabs=0;
end

% %b9 ripples, Appendix B
%8 tl_lambda = tl_ahat*mlambda*nlambda*(1.97-0.44*psihat^.21) ...
%     + ahat*tl_mlambda*nlambda*(1.97-0.44*psihat^.21) ...
%     + ahat*mlambda*tl_nlambda*(1.97-0.44*psihat^.21) ...
%     - ahat*mlambda*nlambda*0.44*.21*psihat^(.21-1)*tl_psihat;
ad_ahat   =ad_ahat   + mlambda*nlambda*(1.97-0.44*psihat^.21)      *ad_lambda;
ad_mlambda=ad_mlambda+ ahat*nlambda*(1.97-0.44*psihat^.21)         *ad_lambda;
ad_nlambda=ad_nlambda+ ahat*mlambda*(1.97-0.44*psihat^.21)         *ad_lambda;
ad_psihat =ad_psihat - ahat*mlambda*nlambda*0.44*.21*psihat^(.21-1)*ad_lambda;
ad_lambda=0;
%7 tl_eta = tl_ahat*meta*neta*(.275-.022*psihat^.42) ...
%          + ahat*tl_meta*neta*(.275-.022*psihat^.42) ...
%          + ahat*meta*tl_neta*(.275-.022*psihat^.42) ...
%          - ahat*meta*neta*.022*.42*psihat^(.42-1)*tl_psihat;
ad_ahat  =ad_ahat  + meta*neta*(.275-.022*psihat^.42)      *ad_eta;
ad_meta  =ad_meta  + ahat*neta*(.275-.022*psihat^.42)      *ad_eta;
ad_neta  =ad_neta  + ahat*meta*(.275-.022*psihat^.42)      *ad_eta;
ad_psihat=ad_psihat- ahat*meta*neta*.022*.42*psihat^(.42-1)*ad_eta;
ad_eta=0;
%6 tl_nlambda=tl_neta;
ad_neta=ad_neta+ad_nlambda;
ad_nlambda=0;
if(psihat<=190)
  %   tl_neta=0;
  ad_neta=0;
elseif(190<psihat&psihat<240)
  %   tl_neta = -.5*sin(pi*(psihat-190)/(240-190)) ...
  %             *(pi*tl_psihat/(240-190));
  ad_psihat=ad_psihat-.5*sin(pi*(psihat-190)/(240-190))*pi/(240-190)*ad_neta;
  ad_neta=0;
else
  %   tl_neta=0;
  ad_neta=0;
end
if(psihatc>psihatt)
  %   tl_psihat=tl_psihatc;
  ad_psihatc=ad_psihatc+ad_psihat;
  ad_psihat=0;
else
  %   tl_psihat=tl_psihatt;
  ad_psihatt=ad_psihatt+ad_psihat;
  ad_psihat=0;
end
%3 tl_psihatt = 1.27^2*2*uhatt*tl_uhatt/((s-1)*g*d50) ...
%     - (1.27*uhatt)^2/((s-1)*g*d50^2)*tl_d50;
ad_uhatt=ad_uhatt+ 1.27^2*2*uhatt/((s-1)*g*d50)  *ad_psihatt;
ad_d50  =ad_d50  - (1.27*uhatt)^2/((s-1)*g*d50^2)*ad_psihatt;
ad_psihatt=0;
%2 tl_psihatc = 1.27^2*2*uhatc*tl_uhatc/((s-1)*g*d50) ...
%     - (1.27*uhatc)^2/((s-1)*g*d50^2)*tl_d50;
ad_uhatc=ad_uhatc+ 1.27^2*2*uhatc/((s-1)*g*d50)  *ad_psihatc;
ad_d50  =ad_d50  - (1.27*uhatc)^2/((s-1)*g*d50^2)*ad_psihatc;
ad_psihatc=0;
if(d50<.22e-3)
  %   tl_meta=0;
  %   tl_mlambda=0;
  ad_meta=0;
  ad_mlambda=0;
elseif(.22e-3<=d50&d50<0.3e-3)
  %   tl_mlambda=.27*tl_d50/(.3e-3-.22e-3);
  ad_d50=ad_d50+.27/(.3e-3-.22e-3)*ad_mlambda;
  ad_mlambda=0;
  %   tl_meta=.45*tl_d50/(.3e-3-.22e-3);
  ad_d50=ad_d50+.45/(.3e-3-.22e-3)*ad_meta;
  ad_meta=0;
else
  %   tl_meta=0;
  %   tl_mlambda=0;
  ad_meta=0;
  ad_mlambda=0;
end

% %b8 wave velocity moments
%3 tl_utildetr=.5*sqrt(2)*tl_uhatt;
ad_uhatt=ad_uhatt+.5*sqrt(2)*ad_utildetr;
ad_utildetr=0;
%2 tl_utildecr=.5*sqrt(2)*tl_uhatc;
ad_uhatc=ad_uhatc+.5*sqrt(2)*ad_utildecr;
ad_utildecr=0;
%1 tl_ahat = tl_uhat.*T/(2*pi) ...
%           + uhat.*tl_T/(2*pi);
ad_uhat=ad_uhat+ T/(2*pi)   *ad_ahat;
ad_T   =ad_T   + uhat/(2*pi)*ad_ahat;
ad_ahat=0;

% %b7 constant parameter mu (eqn A2)
if(d50<=.15e-3)
  %   tl_mu=0;
  ad_mu=0;
elseif(.15e-3<d50&d50<.2e-3)
  %   tl_mu=-5*tl_d50/.05e-3;
  ad_d50=ad_d50-5/.05e-3*ad_mu;
  ad_mu=0;
else
  %   tl_mu=0;
  ad_mu=0;
end

% %b6 critical shields param, Soulsby
%2 tl_theta_cr = .3/(1+1.2*Dstar)^2*(1.2*tl_Dstar) ...
%     - 0.055*exp(-.02*Dstar)*(-.02*tl_Dstar);
ad_Dstar=ad_Dstar+ .3/(1+1.2*Dstar)^2*1.2     *ad_theta_cr;
ad_Dstar=ad_Dstar+ 0.055*exp(-.02*Dstar)*.02  *ad_theta_cr;
ad_theta_cr=0;
%1 tl_Dstar=(g*(s-1)/nu^2)^(1/3)*tl_d50;
ad_d50=ad_d50+(g*(s-1)/nu^2)^(1/3)*ad_Dstar;
ad_Dstar=0;

% %b5 for crest velocities, best I can do without being too fancy and screwing
% % up the AD code is just to consider the TL velocity at the location of the
% % NL model crest/trough
[~,ic]=max(uw);
[~,it]=min(uw);
%2 tl_uhatt = tl_uw(it);
ad_uw(it)=ad_uw(it)+ad_uhatt;
ad_uhatt=0;
%1 tl_uhatc = tl_uw(ic);
ad_uw(ic)=ad_uw(ic)+ad_uhatc;
ad_uhatc=0;

% %b3-4 timing of wave velocity direction change, crest, and trough, based on
% % Ruessink et al 2012.
asinarg=-r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2));
%7 tl_Ttu=tl_Tt-tl_Tc;  % duration of deceleration under trough
ad_Tt=ad_Tt+ad_Ttu;
ad_Tc=ad_Tc-ad_Ttu;
ad_Ttu=0;
%6 tl_Ttu = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,...
%                                       omega,r_r2012,phi_r2012,Ttu);
%5 tl_Tcu = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,...
%                                       omega,r_r2012,phi_r2012,Tcu);
[ad1_omega,ad1_r_r2012,ad1_phi_r2012] = ad_Uwave_ruessink2012_tcrest(ad_Ttu,omega,r_r2012,phi_r2012,Ttu);
[ad2_omega,ad2_r_r2012,ad2_phi_r2012] = ad_Uwave_ruessink2012_tcrest(ad_Tcu,omega,r_r2012,phi_r2012,Tcu);
ad_r_r2012  =ad1_r_r2012  +ad2_r_r2012  ;
ad_phi_r2012=ad1_phi_r2012+ad2_phi_r2012;
ad_omega=ad_omega+ad1_omega+ad2_omega;
%4 tl_T=tl_Tc+tl_Tt;
ad_Tc=ad_Tc+ad_T;
ad_Tt=ad_Tt+ad_T;
ad_T=0;
%3 tl_Tt=-tl_Tc;
ad_Tc=ad_Tc-ad_Tt;
ad_Tt=0;
%2 tl_Tc = 1/sqrt(1-asinarg^2)*tl_asinarg/omega ...
%        - asin(asinarg)/omega^2*tl_omega;
ad_asinarg=ad_asinarg+ 1/sqrt(1-asinarg^2)/omega*ad_Tc;
ad_omega=ad_omega- asin(asinarg)/omega^2*ad_Tc;
ad_Tc=0;
%1 tl_asinarg = -tl_r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)) ...
%     - r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))*tl_phi_r2012 ...
%     + r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)).^2.*( ...
%         .5./sqrt(1-r_r2012^2)*(-2*r_r2012*tl_r_r2012) );
coef=r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)).^2;
ad_r_r2012  =ad_r_r2012  - sin(phi_r2012)/(1+sqrt(1-r_r2012^2))        *ad_asinarg;
ad_phi_r2012=ad_phi_r2012- r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))*ad_asinarg;
ad_r_r2012  =ad_r_r2012  - coef*.5./sqrt(1-r_r2012^2)*2*r_r2012        *ad_asinarg;
ad_asinarg=0;

% %b2 intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% % relevant to van Der A
uw2mean=mean(uw.^2);
%3 tl_uhat=.5./sqrt(2*uw2mean)*2.*tl_uw2mean;
ad_uw2mean=ad_uw2mean+.5./sqrt(2*uw2mean)*2.*ad_uhat;
ad_uhat=0;
%2 tl_uw2mean=mean(2*uw.*tl_uw);
ad_uw = ad_uw + 2*uw*ad_uw2mean/nt;
%1 [tl_uw,tl_r_r2012,tl_phi_r2012]=tl_Uwave_ruessink2012(tl_Hmo,tl_kabs,tl_h,omega*t,uwave_wksp);
[ad1_Hmo,ad1_kabs,ad1_omega,ad1_h]=ad_Uwave_ruessink2012(ad_uw,ad_r_r2012,ad_phi_r2012,omega*t,bkgd.uwave_wksp);
ad_Hmo =ad_Hmo +ad1_Hmo ;
ad_kabs=ad_kabs+ad1_kabs;
ad_omega=ad_omega+ad1_omega;
ad_h   =ad_h   +ad1_h   ;

% %b1 derived params
%2 tl_c=-omega/kabs^2*tl_kabs + tl_omega/kabs;
ad_kabs=ad_kabs- omega/kabs^2*ad_c;
ad_omega = ad_omega + ad_c/kabs;
ad_c=0;
%1 tl_Hmo=tl_Hrms*1.4;
ad_Hrms=ad_Hrms+1.4*ad_Hmo;
ad_Hmo=0;

end  % end of main function
