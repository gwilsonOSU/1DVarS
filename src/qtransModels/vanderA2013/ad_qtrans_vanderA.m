function [ad_d50,ad_d90,ad_h,ad_tanbeta,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_delta,ad_ws,ad_Aw,ad_Sw,ad_Uw,ad_param] = ...
    ad_qtrans_vanderA(ad_qs,bkgd)%,invar)
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
ad_tanbeta    =zeros(nx,1);
ad_Hrms       =zeros(nx,1);
ad_kabs       =zeros(nx,1);
ad_udelta     =zeros(nx,2);
ad_delta      =zeros(nx,1);
ad_Aw         =zeros(nx,1);
ad_Sw         =zeros(nx,1);
ad_Uw         =zeros(nx,1);
ad_omega=0;
ad_tauwRe=0;
ad_streamingEffect=0;
ad_param.n    =0;
ad_param.m    =0;
ad_param.xi   =0;
ad_param.alpha=0;
if(isfield(bkgd(1).param,'Cc'))
  ad_param.Cc   =0;
  ad_param.Cf   =0;
end
ad_param=repmat(ad_param,[nx 1]);

% this wrapper loop serves to handle vector inputs
for i=nx:-1:1
  % tl_qs(i) = ...
  %     tl_qtrans_vanderA(tl_d50,tl_d90,tl_h(i),tl_tanbeta(i),tl_Hrms(i),tl_kabs(i),...
  %                       tl_udelta(i,:),tl_ws,tl_dAw(i),tl_dSw(i),tl_param,bkgd_qtrans(i));
  [ad1_d50,ad1_d90,ad1_h,ad1_tanbeta,ad1_Hrms,ad1_kabs,...
   ad1_omega,ad1_udelta,ad1_delta,ad1_ws,ad1_Aw,ad1_Sw,ad1_Uw,ad1_param] = ...
      ad_qtrans_vanderA_main(ad_qs(i),bkgd(i));%,invar);
  ad_d50(i)=ad_d50(i)+ad1_d50   ;
  ad_d90(i)=ad_d90(i)+ad1_d90   ;
  ad_ws(i) =ad_ws(i) +ad1_ws    ;
  ad_h(i)=ad_h(i)+ad1_h;
  ad_tanbeta(i)=ad_tanbeta(i)+ad1_tanbeta;
  ad_Hrms(i)  =ad_Hrms(i)  +ad1_Hrms  ;
  ad_kabs(i)  =ad_kabs(i)  +ad1_kabs  ;
  ad_omega=ad_omega+ad1_omega;
  ad_udelta(i,:)=ad_udelta(i,:)+ad1_udelta;
  ad_delta(i)=ad_delta(i)+ad1_delta;
  ad_Aw(i)=ad_Aw(i)+ad1_Aw;
  ad_Sw(i)=ad_Sw(i)+ad1_Sw;
  ad_Uw(i)=ad_Uw(i)+ad1_Uw;
  ad_param(i).n    =ad_param(i).n    +ad1_param.n    ;
  ad_param(i).m    =ad_param(i).m    +ad1_param.m    ;
  ad_param(i).xi   =ad_param(i).xi   +ad1_param.xi   ;
  ad_param(i).alpha=ad_param(i).alpha+ad1_param.alpha;
  if(isfield(bkgd(1).param,'Cc'))
    ad_param(i).Cc   =ad_param(i).alpha+ad1_param.Cc   ;
    ad_param(i).Cf   =ad_param(i).alpha+ad1_param.Cf   ;
  end
end

if(~eparam)
  adp2=struct;
  adp2.n    =sum([ad_param.n    ]);
  adp2.m    =sum([ad_param.m    ]);
  adp2.xi   =sum([ad_param.xi   ]);
  adp2.alpha=sum([ad_param.alpha]);
  if(isfield(bkgd(1).param,'Cc'))
    adp2.Cc   =sum([ad_param.Cc   ]);
    adp2.Cf   =sum([ad_param.Cf   ]);
  end
  ad_param=adp2;
end

end  % end of wrapper function, start of main function

function [ad_d50,ad_d90,ad_h,ad_tanbeta,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_delta,ad_ws,ad_Aw,ad_Sw,ad_Uw,ad_param]=ad_qtrans_vanderA_main(ad_qs,bkgd)%,invar)

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
Uw          =bkgd.Uw          ;
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
eps_s        =bkgd.eps_s          ;
uwmo         =bkgd.uwmo           ;
tanbeta      =bkgd.tanbeta        ;
if(isfield(param,'Cc'))  % option-1 for above-WBL transport
  qs2          =bkgd.qs2            ;
  qs3          =bkgd.qs3            ;
else  % option-2 for above-WBL transport
  Lt     =bkgd.Lt      ;
  Lc     =bkgd.Lc      ;
  wfracc =bkgd.wfracc  ;
  wfract =bkgd.wfract  ;
  utot   =bkgd.utot    ;
  Omegas =bkgd.Omegas  ;
  Kdenom1=bkgd.Kdenom1 ;
  Kdenom2=bkgd.Kdenom2 ;
  Kdenom =bkgd.Kdenom  ;
  K      =bkgd.K       ;
end
qsVdA        =bkgd.qsVdA          ;
qsCc         =bkgd.qsCc           ;
qsCf         =bkgd.qsCf           ;

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
ad_Aw   =0;
ad_Sw   =0;
ad_Uw   =0;
ad_param.m    =0;
ad_param.n    =0;
ad_param.alpha=0;
ad_param.xi   =0;
if(isfield(bkgd.param,'Cc'))
  ad_param.Cc   =0;
  ad_param.Cf   =0;
end
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
ad_streamingEffect=0;
ad_tauwRe=0;
ad_argc=0;
ad_argt=0;
ad_argc1=0;
ad_argt1=0;
ad_argc2=0;
ad_argt2=0;
ad_b     =0;
ad_RR    =0;
ad_worb1c=0;
ad_worb1t=0;
ad_worb2c=0;
ad_worb2t=0;
ad_t1ca  =0;
ad_t1ta  =0;
ad_t1cb  =0;
ad_t1tb  =0;
ad_t1c   =0;
ad_t1t   =0;
ad_t2cb  =0;
ad_t2tb  =0;
ad_t2ca  =0;
ad_t2ta  =0;
ad_t2c   =0;
ad_t2t   =0;
ad_worbc =0;
ad_worbt =0;
if(~isfield(param,'streamingType') || param.streamingType=='v')
  ad_fwd=0;
elseif(param.streamingType=='n')
  ad_fws=0;
  ad_r=0;
  ad_f25=0;
  ad_theta25=0;
end
ad_uwmo   =0;
ad_arg_qs2=zeros(1,nt);
ad_qs2    =0;
ad_arg_qs3=zeros(1,nt);
ad_qs3    =0;
ad_qsCc   =0;
ad_qsCf   =0;
ad_qsVdA   =0;
ad_tanbeta=0;
ad_tcr=0;
ad_ttr=0;
ad_tuc=0;
ad_tdc=0;
ad_phidc=0;
ad_phiuc=0;
ad_delta=0;
ad_K=0;
ad_Kdenom=0;
ad_Kdenom1=0;
ad_Kdenom2=0;
ad_arg_Kdenom1=0;
ad_arg_Kdenom2=0;
ad_arg_qsCf=0;
ad_arg_qsCc=0;
ad_utot=zeros(1,nt);
ad_Omegas=0;
ad_wfracc=0;
ad_wfract=0;
ad_Lc=0;
ad_Lt=0;

% % TEST-CODE: override input variable
% if(~strcmp(invar,'qs'))
%   eval(['ad_' invar '=ad_qs;'])
%   ad_qs=0;
% end

% tl_qs = tl_qsVdA + tl_qsCc + tl_qsCf;
ad_qsVdA=ad_qsVdA+ ad_qs;
ad_qsCc =ad_qsCc + ad_qs;
ad_qsCf =ad_qsCf + ad_qs;
ad_qs=0;

% Add suspended load contribution, based on energetics model
term3 = sqrt((s-1)*g*d50^3)/(1-psed);  % NL helper variable
if(isfield(param,'Cc'))  % OPTION-1

  %a8 tl_qsVdA = (tl_qsc + tl_qst)/T*term3 ...
  %       - (qsc + qst)/T^2*term3*tl_T ...
  %       + (qsc + qst)/T*tl_term3;
  ad_qsc  =ad_qsc  + 1/T*term3            *ad_qsVdA;
  ad_qst  =ad_qst  + 1/T*term3            *ad_qsVdA;
  ad_T    =ad_T    - (qsc + qst)/T^2*term3*ad_qsVdA;
  ad_term3=ad_term3+ (qsc + qst)/T        *ad_qsVdA;
  ad_qsVdA=0;
  %a7 tl_qsCf = - tl_qs3*eps_s^2/ws^2*param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
  %           + 2*qs3*eps_s^2/ws^3*tl_ws*param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
  %           - qs3*eps_s^2/ws^2*tl_param.Cf*tanbeta/(g*(s-1)*(1-psed)) ...
  %           - qs3*eps_s^2/ws^2*param.Cf*tl_tanbeta/(g*(s-1)*(1-psed));
  % 
  ad_qs3         =ad_qs3     - eps_s^2/ws^2*param.Cf*tanbeta/(g*(s-1)*(1-psed))      *ad_qsCf;
  ad_ws          =ad_ws      + 2*qs3*eps_s^2/ws^3*param.Cf*tanbeta/(g*(s-1)*(1-psed))*ad_qsCf;
  ad_param.Cf   =ad_param.Cf - qs3*eps_s^2/ws^2*tanbeta/(g*(s-1)*(1-psed))           *ad_qsCf;
  ad_tanbeta     =ad_tanbeta - qs3*eps_s^2/ws^2*param.Cf/(g*(s-1)*(1-psed))          *ad_qsCf;
  ad_qsCf=0;
  %a6 tl_qsCc = + tl_qs2*eps_s/ws*param.Cc     /(g*(s-1)*(1-psed)) ...
  %           - qs2*eps_s/ws^2*param.Cc*tl_ws/(g*(s-1)*(1-psed)) ...
  %           + qs2*eps_s/ws*tl_param.Cc     /(g*(s-1)*(1-psed));
  ad_qs2     =ad_qs2      + eps_s/ws*param.Cc      /(g*(s-1)*(1-psed))*ad_qsCc;
  ad_ws      =ad_ws       - qs2*eps_s/ws^2*param.Cc/(g*(s-1)*(1-psed))*ad_qsCc;
  ad_param.Cc=ad_param.Cc + qs2*eps_s/ws           /(g*(s-1)*(1-psed))*ad_qsCc;
  ad_qsCc=0;
  %a5 tl_qs3 = mean(tl_arg_qs3);
  ad_arg_qs3 = ad_arg_qs3 + 1/nt*ad_qs3*ones(1,nt);  % distribute mean across 1xnt array
  %a4 tl_arg_qs3 = 5*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2).*( uwmo.*tl_uwmo + udelta(1)*tl_udelta(1) + udelta(2)*tl_udelta(2) );
  ad_uwmo     =ad_uwmo     + 5*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2).*uwmo    .*ad_arg_qs3;
  ad_udelta(1)=ad_udelta(1)+ sum(5*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)*udelta(1).*ad_arg_qs3);
  ad_udelta(2)=ad_udelta(2)+ sum(5*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)*udelta(2).*ad_arg_qs3);
  ad_arg_qs3=0;
  %a3 tl_qs2 = mean(tl_arg_qs2);
  ad_arg_qs2 = ad_arg_qs2 + 1/nt*ad_qs2*ones(1,nt);  % distribute mean across 1xnt array
  %a2 tl_arg_qs2 = 3*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(1/2)*udelta(1).*( ...
  %     uwmo.*tl_uwmo + udelta(1)*tl_udelta(1) + udelta(2)*tl_udelta(2) ) ...
  %     + (uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)*tl_udelta(1);
  ad_uwmo     =ad_uwmo     + 3*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(1/2)*udelta(1).*uwmo    .*ad_arg_qs2;
  ad_udelta(1)=ad_udelta(1)+ sum(3*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(1/2)*udelta(1)*udelta(1).*ad_arg_qs2);
  ad_udelta(2)=ad_udelta(2)+ sum(3*(uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(1/2)*udelta(1)*udelta(2).*ad_arg_qs2);
  ad_udelta(1)=ad_udelta(1)+ sum((uwmo.^2 + udelta(1)^2 + udelta(2)^2).^(3/2)                      .*ad_arg_qs2);
  ad_arg_qs2=0;
  %a1 tl_uwmo=tl_uw/1.4;
  ad_uw=ad_uw+1/1.4*ad_uwmo;
  ad_uwmo=0;

else  % OPTION-2

  % NL helper variables
  numst=exp(-deltast/Lt)/0.08;
  numsc=exp(-deltasc/Lc)/0.08;
  denomst=1-Omegat*d50/0.08*deltast/Lt^2*exp(-deltast/Lt);
  denomsc=1-Omegac*d50/0.08*deltasc/Lc^2*exp(-deltasc/Lc);
  dLtdOmegat = numst*d50/denomst;
  dLcdOmegac = numsc*d50/denomsc;
  dLtdd50 = numst*Omegat/denomst;
  dLcdd50 = numsc*Omegac/denomsc;
  arg_Kdenom1=utot.^3;
  arg_Kdenom2 = utot.^4;
  arg_qsCc = utot.^3*udelta(1);
  arg_qsCf = utot.^5;

  %b18 tl_qsCf = ...
  %     - 1/(g*(s-1)*(1-psed))*eps_s^2*tl_K/ws^2*mean(arg_qsCf)*tanbeta ...
  %     + 2/(g*(s-1)*(1-psed))*eps_s^2*K/ws^3*mean(arg_qsCf)*tanbeta*tl_ws ...
  %     - 1/(g*(s-1)*(1-psed))*eps_s^2*K/ws^2*mean(tl_arg_qsCf)*tanbeta ...
  %     - 1/(g*(s-1)*(1-psed))*eps_s^2*K/ws^2*mean(arg_qsCf)*tl_tanbeta;
  ad_K       =ad_K       - 1/(g*(s-1)*(1-psed))*eps_s^2/ws^2*mean(arg_qsCf)*tanbeta   *ad_qsCf;
  ad_ws      =ad_ws      + 2/(g*(s-1)*(1-psed))*eps_s^2*K/ws^3*mean(arg_qsCf)*tanbeta *ad_qsCf;
  ad_arg_qsCf=ad_arg_qsCf- 1/(g*(s-1)*(1-psed))*eps_s^2*K/ws^2*tanbeta*1/nt*ones(1,nt)*ad_qsCf;  % distribute mean across 1xnt array
  ad_tanbeta =ad_tanbeta - 1/(g*(s-1)*(1-psed))*eps_s^2*K/ws^2*mean(arg_qsCf)         *ad_qsCf;
  ad_qsCf=0;
  %b17 tl_arg_qsCf = 5*utot.^4.*tl_utot;
  ad_utot=ad_utot+ 5*utot.^4.*ad_arg_qsCf;
  %b16 tl_qsCc = ...
  %     + 1/(g*(s-1)*(1-psed))*eps_s*tl_K/ws*mean(arg_qsCc) ...
  %     - 1/(g*(s-1)*(1-psed))*eps_s*K/ws^2*mean(arg_qsCc)*tl_ws ...
  %     + 1/(g*(s-1)*(1-psed))*eps_s*K/ws*mean(tl_arg_qsCc);
  ad_K       =ad_K       + 1/(g*(s-1)*(1-psed))*eps_s/ws*mean(arg_qsCc)                    *ad_qsCc;
  ad_ws      =ad_ws      - 1/(g*(s-1)*(1-psed))*eps_s*K/ws^2*mean(arg_qsCc)                *ad_qsCc;
  ad_arg_qsCc=ad_arg_qsCc+ 1/(g*(s-1)*(1-psed))*eps_s*K/ws                 *1/nt*ones(1,nt)*ad_qsCc;  % distribute mean across 1xnt array
  ad_qsCc=0;
  %b15 tl_arg_qsCc = + 3*utot.^2*udelta(1).*tl_utot ...
  %     + utot.^3*tl_udelta(1);
  ad_utot     =ad_utot     + 3*utot.^2*udelta(1).*ad_arg_qsCc;
  ad_udelta(1)=ad_udelta(1)+ sum(utot.^3        .*ad_arg_qsCc);
  ad_arg_qsCc=0;
  %b14 tl_K = ...
  %     + (s-1)*g*tl_d50*Omegas/(Kdenom1+Kdenom2) ...
  %     + (s-1)*g*d50*tl_Omegas/(Kdenom1+Kdenom2) ...
  %     - (s-1)*g*d50*Omegas/(Kdenom1+Kdenom2)^2*( tl_Kdenom1 + tl_Kdenom2 );
  ad_d50    =ad_d50    + (s-1)*g*Omegas/(Kdenom1+Kdenom2)      *ad_K;
  ad_Omegas =ad_Omegas + (s-1)*g*d50/(Kdenom1+Kdenom2)         *ad_K;
  ad_Kdenom1=ad_Kdenom1- (s-1)*g*d50*Omegas/(Kdenom1+Kdenom2)^2*ad_K;
  ad_Kdenom2=ad_Kdenom2- (s-1)*g*d50*Omegas/(Kdenom1+Kdenom2)^2*ad_K;
  ad_K=0;
  %b13 tl_Kdenom = tl_Kdenom1 + tl_Kdenom2;
  ad_Kdenom1=ad_Kdenom1+ ad_Kdenom;
  ad_Kdenom2=ad_Kdenom2+ ad_Kdenom;
  ad_Kdenom=0;
  %b12 tl_Kdenom2 = ...
  %     + 2*eps_s^2/ws^3*mean(arg_Kdenom2)*tanbeta*tl_ws ...
  %     - eps_s^2/ws^2*mean(tl_arg_Kdenom2)*tanbeta ...
  %     - eps_s^2/ws^2*mean(arg_Kdenom2)*tl_tanbeta;
  ad_ws         =ad_ws         + 2*eps_s^2/ws^3*mean(arg_Kdenom2)*tanbeta*ad_Kdenom2;
  ad_arg_Kdenom2=ad_arg_Kdenom2- eps_s^2/ws^2*tanbeta*1/nt*ones(1,nt)    *ad_Kdenom2;  % distribute mean across 1xnt array
  ad_tanbeta    =ad_tanbeta    - eps_s^2/ws^2*mean(arg_Kdenom2)          *ad_Kdenom2;
  ad_Kdenom2=0;
  %b11 tl_arg_Kdenom2 = 4*utot.^3.*tl_utot;
  ad_utot=ad_utot+ 4*utot.^3.*ad_arg_Kdenom2;
  ad_arg_Kdenom2=0;
  %b10 tl_Kdenom1 = ...
  %     - eps_s/ws^2*mean(arg_Kdenom1)*tl_ws ...
  %     + eps_s/ws*mean(tl_arg_Kdenom1);
  ad_ws         =ad_ws         - eps_s/ws^2*mean(arg_Kdenom1)                *ad_Kdenom1;
  ad_arg_Kdenom1=ad_arg_Kdenom1+ eps_s/ws                    *1/nt*ones(1,nt)*ad_Kdenom1;  % distribute mean across 1xnt array
  ad_Kdenom1=0;
  %b9 tl_arg_Kdenom1 = 3*utot.^2.*tl_utot;
  ad_utot=ad_utot+ 3*utot.^2.*ad_arg_Kdenom1;
  ad_arg_Kdenom1=0;
  %b8 tl_utot = 1./sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).*( ...
  %     + uwmo     .*tl_uwmo      ...
  %     + udelta(1)*tl_udelta(1) ...
  %     + udelta(2)*tl_udelta(2) );
  ad_uwmo     =ad_uwmo     + 1./sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).*uwmo     .*ad_utot;
  ad_udelta(1)=ad_udelta(1)+ sum(1./sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).*udelta(1).*ad_utot);
  ad_udelta(2)=ad_udelta(2)+ sum(1./sqrt(uwmo.^2+udelta(1)^2+udelta(2)^2).*udelta(2).*ad_utot);
  %b7 tl_Omegas = ...
  %     - tl_wfracc*Omegac*Tc/T ...
  %     + (1-wfracc)*tl_Omegac*Tc/T ...
  %     + (1-wfracc)*Omegac*tl_Tc/T ...
  %     - (1-wfracc)*Omegac*Tc/T^2*tl_T ...
  %     - tl_wfract*Omegat*Tt/T ...
  %     + (1-wfract)*tl_Omegat*Tt/T ...
  %     + (1-wfract)*Omegat*tl_Tt/T ...
  %     - (1-wfract)*Omegat*Tt/T^2*tl_T;
  ad_wfracc=ad_wfracc- Omegac*Tc/T             *ad_Omegas;
  ad_Omegac=ad_Omegac+ (1-wfracc)*Tc/T         *ad_Omegas;
  ad_Tc    =ad_Tc    + (1-wfracc)*Omegac/T     *ad_Omegas;
  ad_T     =ad_T     - (1-wfracc)*Omegac*Tc/T^2*ad_Omegas;
  ad_wfract=ad_wfract- Omegat*Tt/T             *ad_Omegas;
  ad_Omegat=ad_Omegat+ (1-wfract)*Tt/T         *ad_Omegas;
  ad_Tt    =ad_Tt    + (1-wfract)*Omegat/T     *ad_Omegas;
  ad_T     =ad_T     - (1-wfract)*Omegat*Tt/T^2*ad_Omegas;
  ad_Omegas=0;
  %b6 tl_uwmo = + 1/1.4*tl_uw;
  ad_uw=ad_uw+ 1/1.4*ad_uwmo;
  ad_uwmo=0;
  %b5 tl_qsVdA = (tl_qsc*wfracc + qsc*tl_wfracc + tl_qst*wfract + qst*tl_wfract)/T*term3 ...
  %       - (qsc*wfracc + qst*wfract)/T^2*term3*tl_T ...
  %       + (qsc*wfracc + qst*wfract)/T*tl_term3;
  ad_qsc   =ad_qsc    + 1/T*term3*wfracc                   *ad_qsVdA;
  ad_wfracc=ad_wfracc + 1/T*term3*qsc                      *ad_qsVdA;
  ad_qst   =ad_qst    + 1/T*term3*wfract                   *ad_qsVdA;
  ad_wfract=ad_wfract + 1/T*term3*qst                      *ad_qsVdA;
  ad_T     =ad_T      - (qsc*wfracc + qst*wfract)/T^2*term3*ad_qsVdA;
  ad_term3 =ad_term3  + (qsc*wfracc + qst*wfract)/T        *ad_qsVdA;
  ad_qsVdA=0;
  %b4 tl_wfract = - delta/Lt^2*exp(-delta/Lt)*tl_Lt ...
  %     + 1/Lt*exp(-delta/Lt)*tl_delta;
  ad_Lt   =ad_Lt   - delta/Lt^2*exp(-delta/Lt)*ad_wfract;
  ad_delta=ad_delta+ 1/Lt*exp(-delta/Lt)      *ad_wfract;
  ad_wfract=0;
  %b3 tl_wfracc = - delta/Lc^2*exp(-delta/Lc)*tl_Lc ...
  %     + 1/Lc*exp(-delta/Lc)*tl_delta;
  ad_Lc   =ad_Lc   - delta/Lc^2*exp(-delta/Lc)*ad_wfracc;
  ad_delta=ad_delta+ 1/Lc*exp(-delta/Lc)      *ad_wfracc;
  ad_wfracc=0;
  %b2 tl_Lc = + dLcdOmegac*tl_Omegac ...
  %         + dLcdd50*tl_d50;
  ad_Omegac=ad_Omegac+ dLcdOmegac*ad_Lc;
  ad_d50   =ad_d50   + dLcdd50   *ad_Lc;
  ad_Lc=0;
  %b1 tl_Lt = + dLtdOmegat*tl_Omegat ...
  %         + dLtdd50*tl_d50;
  ad_Omegat=ad_Omegat+ dLtdOmegat*ad_Lt;
  ad_d50   =ad_d50   + dLtdd50   *ad_Lt;
  ad_Lt=0;

end

%b14 transport, eqn 1
absthetac=abs(thetac);
absthetat=abs(thetat);
term3 = sqrt((s-1)*g*d50^3)/(1-psed);
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

%18 tl_wst = tl_ws + tl_worbt;
ad_ws   =ad_ws   + ad_wst;
ad_worbt=ad_worbt+ ad_wst;
ad_wst=0;
%17 tl_wsc = tl_ws - tl_worbc;
ad_ws   =ad_ws   + ad_wsc;
ad_worbc=ad_worbc- ad_wsc;
ad_wsc=0;

% Use stokes 2nd order theory to get vertical fluid velocities
%20 tl_worbt = .125*tl_worb1t*sqrt(t1t) ...
%     + .5*.125*worb1t/sqrt(t1t)*tl_t1t ...
%     + tl_worb2t*sin(t2t) ...
%     + worb2t*cos(t2t)*tl_t2t;
ad_worb1t=ad_worb1t+ .125*sqrt(t1t)          *ad_worbt;
ad_t1t   =ad_t1t   + .5*.125*worb1t/sqrt(t1t)*ad_worbt;
ad_worb2t=ad_worb2t+ sin(t2t)                *ad_worbt;
ad_t2t   =ad_t2t   + worb2t*cos(t2t)         *ad_worbt;
ad_worbt=0;
%19 tl_worbc = .125*tl_worb1c*sqrt(t1c) ...
%     + .5*.125*worb1c/sqrt(t1c)*tl_t1c ...
%     + tl_worb2c*sin(t2c) ...
%     + worb2c*cos(t2c)*tl_t2c;
ad_worb1c=ad_worb1c+ .125*sqrt(t1c)          *ad_worbc;
ad_t1c   =ad_t1c   + .5*.125*worb1c/sqrt(t1c)*ad_worbc;
ad_worb2c=ad_worb2c+ sin(t2c)                *ad_worbc;
ad_t2c   =ad_t2c   + worb2c*cos(t2c)         *ad_worbc;
ad_worbc=0;
%18 tl_t2t = -2/sqrt(1-t2ta^2)*tl_t2ta;
ad_t2ta=ad_t2ta- 2/sqrt(1-t2ta^2)*ad_t2t;
ad_t2t=0;
%17 tl_t2c = -2/sqrt(1-t2ca^2)*tl_t2ca;
ad_t2ca=ad_t2ca- 2/sqrt(1-t2ca^2)*ad_t2c;
ad_t2c=0;
%16 tl_t2ta = .125*tl_t2tb/worb2t ...
%           - .125*t2tb/worb2t^2*tl_worb2t;
ad_t2tb  =ad_t2tb  + .125/worb2t       *ad_t2ta;
ad_worb2t=ad_worb2t- .125*t2tb/worb2t^2*ad_t2ta;
ad_t2ta=0;
%15 tl_t2ca = .125*tl_t2cb/worb2c ...
%           - .125*t2cb/worb2c^2*tl_worb2c;
ad_t2cb  =ad_t2cb  + .125/worb2c       *ad_t2ca;
ad_worb2c=ad_worb2c- .125*t2cb/worb2c^2*ad_t2ca;
ad_t2ca=0;
%14 tl_t2tb = -tl_worb1t + 1/sqrt( worb1t^2 + 32*worb2t^2 )*( worb1t*tl_worb1t + 32*worb2t*tl_worb2t );
ad_worb1t=ad_worb1t- 1                                         *ad_t2tb;
ad_worb1t=ad_worb1t+ 1/sqrt( worb1t^2 + 32*worb2t^2 )*worb1t   *ad_t2tb;
ad_worb2t=ad_worb2t+ 1/sqrt( worb1t^2 + 32*worb2t^2 )*32*worb2t*ad_t2tb;
ad_t2tb=0;
%13 tl_t2cb = -tl_worb1c + 1/sqrt( worb1c^2 + 32*worb2c^2 )*( worb1c*tl_worb1c + 32*worb2c*tl_worb2c );
ad_worb1c=ad_worb1c- 1                                         *ad_t2cb;
ad_worb1c=ad_worb1c+ 1/sqrt( worb1c^2 + 32*worb2c^2 )*worb1c   *ad_t2cb;
ad_worb2c=ad_worb2c+ 1/sqrt( worb1c^2 + 32*worb2c^2 )*32*worb2c*ad_t2cb;
ad_t2cb=0;
%12 tl_t1t = 2*t1tb/worb2t^2*tl_t1tb ...
%          - 2*t1tb^2/worb2t^3*tl_worb2t;
ad_t1tb  =ad_t1tb  + 2*t1tb/worb2t^2  *ad_t1t;
ad_worb2t=ad_worb2t- 2*t1tb^2/worb2t^3*ad_t1t;
ad_t1t=0;
%11 tl_t1c = 2*t1cb/worb2c^2*tl_t1cb ...
%          - 2*t1cb^2/worb2c^3*tl_worb2c;
ad_t1cb  =ad_t1cb  + 2*t1cb/worb2c^2  *ad_t1c;
ad_worb2c=ad_worb2c- 2*t1cb^2/worb2c^3*ad_t1c;
ad_t1c=0;
%10 tl_t1tb = tl_worb1t - .5/sqrt(t1ta)*tl_t1ta;
ad_worb1t=ad_worb1t+ 1            *ad_t1tb;
ad_t1ta  =ad_t1ta  - .5/sqrt(t1ta)*ad_t1tb;
ad_t1tb=0;
%9 tl_t1cb = tl_worb1c - .5/sqrt(t1ca)*tl_t1ca;
ad_worb1c=ad_worb1c+ 1            *ad_t1cb;
ad_t1ca  =ad_t1ca  - .5/sqrt(t1ca)*ad_t1cb;
ad_t1cb=0;
%8 tl_t1ta = 2*worb1t*tl_worb1t + 32*2*worb2t*tl_worb2t;
ad_worb1t=ad_worb1t+ 2*worb1t   *ad_t1ta;
ad_worb2t=ad_worb2t+ 32*2*worb2t*ad_t1ta;
ad_t1ta=0;
%7 tl_t1ca = 2*worb1c*tl_worb1c + 32*2*worb2c*tl_worb2c;
ad_worb1c=ad_worb1c+ 2*worb1c   *ad_t1ca;
ad_worb2c=ad_worb2c+ 32*2*worb2c*ad_t1ca;
ad_t1ca=0;
%6 tl_worb2t = 2*tl_worb1t*(2*RR-1) ...
%     + 2*worb1t*2*tl_RR;
ad_worb1t=ad_worb1t+ 2*(2*RR-1)*ad_worb2t;
ad_RR    =ad_RR    + 2*worb1t*2*ad_worb2t;
ad_worb2t=0;
%5 tl_worb2c = 2*tl_worb1c*(2*RR-1) ...
%     + 2*worb1c*2*tl_RR;
ad_worb1c=ad_worb1c+ 2*(2*RR-1)*ad_worb2c;
ad_RR    =ad_RR    + 2*worb1c*2*ad_worb2c;
ad_worb2c=0;
%4 tl_worb1t = + pi*tl_Hmo*etawt/(T*h) ...
%     + pi*Hmo*tl_etawt/(T*h) ...
%     - pi*Hmo*etawt/(T^2*h)*tl_T ...
%     - pi*Hmo*etawt/(T*h^2)*tl_h;
ad_Hmo  =ad_Hmo  + pi*etawt/(T*h)      *ad_worb1t;
ad_etawt=ad_etawt+ pi*Hmo/(T*h)        *ad_worb1t;
ad_T    =ad_T    - pi*Hmo*etawt/(T^2*h)*ad_worb1t;
ad_h    =ad_h    - pi*Hmo*etawt/(T*h^2)*ad_worb1t;
ad_worb1t=0;
%3 tl_worb1c = + pi*tl_Hmo*etawc/(T*h) ...
%     + pi*Hmo*tl_etawc/(T*h) ...
%     - pi*Hmo*etawc/(T^2*h)*tl_T ...
%     - pi*Hmo*etawc/(T*h^2)*tl_h;
ad_Hmo  =ad_Hmo  + pi*etawc/(T*h)      *ad_worb1c;
ad_etawc=ad_etawc+ pi*Hmo/(T*h)        *ad_worb1c;
ad_T    =ad_T    - pi*Hmo*etawc/(T^2*h)*ad_worb1c;
ad_h    =ad_h    - pi*Hmo*etawc/(T*h^2)*ad_worb1c;
ad_worb1c=0;
%2 tl_RR = -0.5*tl_b*sin(phi_r2012) ...
%         -0.5*b*cos(phi_r2012)*tl_phi_r2012;
ad_b        =ad_b        - 0.5*sin(phi_r2012)  *ad_RR;
ad_phi_r2012=ad_phi_r2012- 0.5*b*cos(phi_r2012)*ad_RR;
ad_RR=0;
%1 tl_b = -1/r_r2012^2*(1-sqrt(1-r_r2012^2))*tl_r_r2012 ...
%        + 1/r_r2012/sqrt(1-r_r2012^2)*r_r2012*tl_r_r2012;
ad_r_r2012=ad_r_r2012- 1/r_r2012^2*(1-sqrt(1-r_r2012^2))  *ad_b;
ad_r_r2012=ad_r_r2012+ 1/r_r2012/sqrt(1-r_r2012^2)*r_r2012*ad_b;
ad_b=0;

% %b12 sheet flow layer thickness, Appendix C
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

% boundary layer streaming, either with VDA13 model or Nielsen
alphaw = 4/(3*pi);
if(~isfield(param,'streamingType') || param.streamingType=='v')
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
elseif(param.streamingType=='n')
  %6 tl_tauwRe = tl_fws*alphaw*uhat^3/2/c ...
  %     + 3*fws*alphaw*uhat^2/2/c*tl_uhat ...
  %     - fws*alphaw*uhat^3/2/c^2*tl_c;
  ad_fws   =ad_fws   + alphaw*uhat^3/2/c      *ad_tauwRe;
  ad_uhat  =ad_uhat  + 3*fws*alphaw*uhat^2/2/c*ad_tauwRe;
  ad_c     =ad_c     - fws*alphaw*uhat^3/2/c^2*ad_tauwRe;
  ad_tauwRe=0;
  %5 tl_fws = exp(5.5*(r/ahat)^.2-6.3)*2*5.5*(r/ahat)*( tl_r/ahat - r/ahat.^2.*tl_ahat );
  ad_r   =ad_r   + exp(5.5*(r/ahat)^.2-6.3)*2*5.5*(r/ahat)/ahat      *ad_fws;
  ad_ahat=ad_ahat- exp(5.5*(r/ahat)^.2-6.3)*2*5.5*(r/ahat)*r/ahat.^2.*ad_fws;
  ad_fws=0;
  if(lambda>0)
    %4   tl_r = tl_r + 2*8*eta/lambda*tl_eta - 8*eta^2/lambda^2*tl_lambda;
    ad_eta   =ad_eta   + 2*8*eta/lambda  *ad_r;
    ad_lambda=ad_lambda- 8*eta^2/lambda^2*ad_r;
    % do not clear ad_r, statement was an increment not assignment
  end
  if(theta25-0.05<=0)
    %3 tl_r=0;
    ad_r=0;
  else
    %3 tl_r = 170*sqrt(theta25-0.05)*tl_d50 ...
    %      + .5*170/sqrt(theta25-0.05)*tl_theta25;
    ad_d50    =ad_d50    + 170*sqrt(theta25-0.05)   *ad_r;
    ad_theta25=ad_theta25+ .5*170/sqrt(theta25-0.05)*ad_r;
    ad_r=0;
  end
  %2 tl_theta25 = 0.5*tl_f25*(ahat*omega)^2/((s-1)*g*d50) ...
  %     + 2*0.5*f25*ahat*(omega)^2/((s-1)*g*d50)*tl_ahat ...
  %     + 2*0.5*f25*(ahat)^2*omega/((s-1)*g*d50)*tl_omega ...
  %     - 0.5*f25*(ahat*omega)^2/((s-1)*g*d50^2)*tl_d50;
  ad_f25  =ad_f25  + 0.5*(ahat*omega)^2/((s-1)*g*d50)      *ad_theta25;
  ad_ahat =ad_ahat + 2*0.5*f25*ahat*(omega)^2/((s-1)*g*d50)*ad_theta25;
  ad_omega=ad_omega+ 2*0.5*f25*(ahat)^2*omega/((s-1)*g*d50)*ad_theta25;
  ad_d50  =ad_d50  - 0.5*f25*(ahat*omega)^2/((s-1)*g*d50^2)*ad_theta25;
  ad_theta25=0;
  %1 tl_f25 = exp(5.5*(2.5*d50/ahat)^2-6.3)*2*5.5*(2.5*d50/ahat)*( ...
  %     2.5*tl_d50/ahat - 2.5*d50/ahat^2*tl_ahat );
  ad_d50 =ad_d50 + exp(5.5*(2.5*d50/ahat)^2-6.3)*2*5.5*(2.5*d50/ahat)*2.5/ahat      *ad_f25;
  ad_ahat=ad_ahat- exp(5.5*(2.5*d50/ahat)^2-6.3)*2*5.5*(2.5*d50/ahat)*2.5*d50/ahat^2*ad_f25;
  ad_f25=0;
elseif(param.streamingType=='0')
  % tl_tauwRe=0;
  ad_tauwRe=0;
else
  error('must provide param.streamingType as either ''n'' or ''v''')
end

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

%5b_orig tl_utrvec = tl_utildetr*[-1 0] + tl_udelta;
%5b_rewrite1 tl_utrvec(1) = -tl_utildetr + tl_udelta(1);
ad_utildetr =ad_utildetr - ad_utrvec(1);
ad_udelta(1)=ad_udelta(1)+ ad_utrvec(1);
ad_utrvec(1)=0;
%5b_rewrite2 tl_utrvec(2) = + tl_udelta(2);
ad_udelta(2)=ad_udelta(2)+ad_utrvec(2);
ad_utrvec(2)=0;
%5a_orig tl_ucrvec = tl_utildecr*[+1 0] + tl_udelta;
%5a_rewrite1 tl_ucrvec(1) = + tl_utildecr + tl_udelta(1);
ad_utildecr =ad_utildecr + ad_ucrvec(1);
ad_udelta(1)=ad_udelta(1)+ ad_ucrvec(1);
ad_ucrvec(1)=0;
%5a_rewrite2 tl_ucrvec(2) = + tl_udelta(2);
ad_udelta(2)=ad_udelta(2)+ad_ucrvec(2);
ad_ucrvec(2)=0;

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
%2 tl_theta_cr = -.3/(1+1.2*Dstar)^2*(1.2*tl_Dstar) ...
%     - 0.055*exp(-.02*Dstar)*(-.02*tl_Dstar);
ad_Dstar=ad_Dstar- .3/(1+1.2*Dstar)^2*1.2     *ad_theta_cr;
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
%2 tl_uhatt = -tl_uw(it);
ad_uw(it)=ad_uw(it)-ad_uhatt;
ad_uhatt=0;
%1 tl_uhatc = tl_uw(ic);
ad_uw(ic)=ad_uw(ic)+ad_uhatc;
ad_uhatc=0;

% timing of wave velocity direction change, crest, and trough, based on
% Ruessink et al 2012.
asinarg=-r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2));
%12 tl_Tcu=tl_tcr-tl_tuc;
ad_tcr=ad_tcr+ad_Tcu;
ad_tuc=ad_tuc-ad_Tcu;
ad_Tcu=0;
%11 tl_Ttu=tl_ttr-tl_tdc;
ad_ttr=ad_ttr+ad_Ttu;
ad_tdc=ad_tdc-ad_Ttu;
ad_Ttu=0;
%10 tl_Tt=tl_T-tl_Tc;
ad_T =ad_T +ad_Tt;
ad_Tc=ad_Tc-ad_Tt;
ad_Tt=0;
%9 tl_Tc = tl_tdc-tl_tuc;
ad_tdc=ad_tdc+ad_Tc;
ad_tuc=ad_tuc-ad_Tc;
ad_Tc=0;
%8 tl_T=-2*pi/omega^2*tl_omega;
ad_omega=ad_omega-2*pi/omega^2*ad_T;
ad_T=0;
%7 tl_ttr = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,omega,r_r2012,phi_r2012,t(itu_guess));
[ad1_omega,ad1_r_r2012,ad1_phi_r2012] = ad_Uwave_ruessink2012_tcrest(ad_ttr,omega,r_r2012,phi_r2012,t(itu_guess));
ad_omega    =ad_omega    +ad1_omega    ;
ad_r_r2012  =ad_r_r2012  +ad1_r_r2012  ;
ad_phi_r2012=ad_phi_r2012+ad1_phi_r2012;
%6 tl_tcr = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r_r2012,tl_phi_r2012,omega,r_r2012,phi_r2012,t(icu_guess));
[ad1_omega,ad1_r_r2012,ad1_phi_r2012] = ad_Uwave_ruessink2012_tcrest(ad_tcr,omega,r_r2012,phi_r2012,t(icu_guess));
ad_omega    =ad_omega    +ad1_omega    ;
ad_r_r2012  =ad_r_r2012  +ad1_r_r2012  ;
ad_phi_r2012=ad_phi_r2012+ad1_phi_r2012;
%5 tl_tdc = tl_phidc/omega - phidc/omega^2*tl_omega;
ad_phidc=ad_phidc+ 1/omega      *ad_tdc;
ad_omega=ad_omega- phidc/omega^2*ad_tdc;
ad_tdc=0;
%4 tl_tuc = tl_phiuc/omega - phiuc/omega^2*tl_omega;
ad_phiuc=ad_phiuc+1/omega       *ad_tuc;
ad_omega=ad_omega- phiuc/omega^2*ad_tuc;
ad_tuc=0;
%3 tl_phidc=-tl_phiuc;
ad_phiuc=ad_phiuc-ad_phidc;
ad_phidc=0;
%2 tl_phiuc = 1./sqrt(1-asinarg^2)/omega*tl_asinarg - asinarg/omega^2*tl_omega;
ad_asinarg=ad_asinarg+ 1./sqrt(1-asinarg^2)/omega*ad_phiuc;
ad_omega  =ad_omega  - asinarg/omega^2           *ad_phiuc;
ad_phiuc=0;
%1 tl_asinarg = -tl_r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)) ...
%     - r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))*tl_phi_r2012 ...
%     + r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2)).^2.*( ...
%         .5./sqrt(1-r_r2012^2)*(-2*r_r2012*tl_r_r2012) );
ad_r_r2012  =ad_r_r2012  - sin(phi_r2012)/(1+sqrt(1-r_r2012^2))                                             *ad_asinarg;
ad_phi_r2012=ad_phi_r2012- r_r2012*cos(phi_r2012)/(1+sqrt(1-r_r2012^2))                                     *ad_asinarg;
ad_r_r2012  =ad_r_r2012  + r_r2012*sin(phi_r2012)/(1+sqrt(1-r_r2012^2))^2*.5./sqrt(1-r_r2012^2)*(-2)*r_r2012*ad_asinarg;
ad_asinarg=0;

% %b2 intra-wave velocities using Ruessink et al. (2012).  Then extract stats
% % relevant to van Der A
uw2mean=mean(uw.^2);
%3 tl_uhat=.5./sqrt(2*uw2mean)*2.*tl_uw2mean;
ad_uw2mean=ad_uw2mean+.5./sqrt(2*uw2mean)*2.*ad_uhat;
ad_uhat=0;
%2 tl_uw2mean=mean(2*uw.*tl_uw);
ad_uw = ad_uw + 2*uw*ad_uw2mean/nt;
%1 [tl_uw,tl_r_r2012,tl_phi_r2012]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,omega*t,uwave_wksp);
[ad1_Aw,ad1_Sw,ad1_Uw]=ad_Uwave_ruessink2012(ad_uw,ad_r_r2012,ad_phi_r2012,omega*t,uwave_wksp);
ad_Aw=ad_Aw+ad1_Aw;
ad_Sw=ad_Sw+ad1_Sw;
ad_Uw=ad_Uw+ad1_Uw;

% note, van der A specifies to use significant orbital velocity amplitude,
% while Ruessink et al. (2012) uses rms.  Uwave_ruessink2012() follows the
% Ruessink et al. (2012) convention, so I need to revert to van der A's
% convention here as a special case.
% tl_Uw=1.4*tl_Uw;
ad_Uw=1.4*ad_Uw;

% %b1 derived params
%2 tl_c=-omega/kabs^2*tl_kabs + tl_omega/kabs;
ad_kabs=ad_kabs- omega/kabs^2*ad_c;
ad_omega = ad_omega + ad_c/kabs;
ad_c=0;
%1 tl_Hmo=tl_Hrms*1.4;
ad_Hrms=ad_Hrms+1.4*ad_Hmo;
ad_Hmo=0;

end  % end of main function
