%
% symmetry test
%
clear
addpath(genpath('../../..'));

param.streamingType='v';  % choose either 'n' or 'v'.  NOTE, should test both

% TEST-CODE: choose a test variable.  These are in the order they appear in
% TL model.  This is a very effective debugging technique for the AD model.
% Requires having an auxiliary input for the TL and AD codes that override
% the default output variable.  See {tl,ad}_qtrans_vanderA.m 'TEST-CODE'
% blocks (commented out for production code)
% inoutvar='d50';
% inoutvar='d90';
% inoutvar='h';
% inoutvar='Hrms';
% inoutvar='kabs';
% inoutvar='omega';
% inoutvar='udelta';  % can't do, vector
% inoutvar='ws';
% inoutvar='asinarg';
% inoutvar='Tc';
% inoutvar='Tcu';
% inoutvar='Ttu';
% inoutvar='uhatc';
% inoutvar='uhatt';
% inoutvar='theta_cr';
% inoutvar='ahat';
% inoutvar='psihatc';
% inoutvar='psihatt';
% inoutvar='neta';  % nan
% inoutvar='eta';  % nan
% inoutvar='lambda';  % nan
% inoutvar='theta_av';
% inoutvar='ksd';
% inoutvar='ksw';
% inoutvar='fd';
% inoutvar='fw';
% inoutvar='alpha';
% inoutvar='argc2';
% inoutvar='argc1';
% inoutvar='argc';
% inoutvar='argt2';
% inoutvar='argt1';
% inoutvar='argt';
% inoutvar='fwdc';
% inoutvar='fwdt';
% inoutvar='ucrvec';  % can't do, vector
% inoutvar='utrvec';  % can't do, vector
% inoutvar='ucrabs';
% inoutvar='utrabs';
% inoutvar='thetac';
% inoutvar='thetat';
% inoutvar='fwd';
% inoutvar='Hmo';
% inoutvar='c';
% inoutvar='uw2mean';
% inoutvar='uhat';
% inoutvar='asinarg ';
% inoutvar='T';  % nan
% inoutvar='Tt';
% inoutvar='Dstar';
% inoutvar='mu';
% inoutvar='ahat ';
% inoutvar='utildecr';
% inoutvar='utildetr';
% inoutvar='meta';  %nan
% inoutvar='mlambda'; %nan
% inoutvar='psihatc ';
% inoutvar='psihatt ';
% inoutvar='psihat';
% inoutvar='neta';
% inoutvar='nlambda';
% inoutvar='eta ';
% inoutvar='lambda ';
% inoutvar='udabs';
% inoutvar='fwc';
% inoutvar='fwt';
% inoutvar='etawc';
% inoutvar='etawt';
% inoutvar='phi_r2012';  %nan
% inoutvar='b ';  % nan
% inoutvar='RR ';  % nan
% inoutvar='worb1c ';  % ok
% inoutvar='worb1t '; %ok
% inoutvar='worb2c ';  % ok
% inoutvar='worb2t ';  % ok
% inoutvar='t1ca';
% inoutvar='asinarg '; %ok
% inoutvar='tauwRe';  % ok
% inoutvar='streamingEffect';
% inoutvar='thetacx';
% inoutvar='thetatx';
% inoutvar='thetahatc';
% inoutvar='thetahatt';
% inoutvar='deltasc';
% inoutvar='deltast';
% inoutvar='deltasc';
% inoutvar='deltast';
% inoutvar='b';  % nan
% inoutvar='RR'; % nan
% inoutvar='worb1c';
% inoutvar='worb1t';
% inoutvar='worb2c';  % NO
% inoutvar='worb2t';  % NO
% inoutvar='t1ca';    % NO
% inoutvar='t1c';     % NO
% inoutvar='t2ca';
% inoutvar='t2c';
% inoutvar='worbc';
% inoutvar='t1ta';
% inoutvar='t1t';
% inoutvar='t2ta';
% inoutvar='t2t';
% inoutvar='worbt';
% inoutvar='wsc';
% inoutvar='wst';
% inoutvar='Pc';
% inoutvar='Pt';
% inoutvar='absthetac';
% inoutvar='absthetat';
% inoutvar='Omegac';
% inoutvar='Omegat';
% inoutvar='Omegacc';
% inoutvar='Omegatt';
% inoutvar='Omegact';  % nan
% inoutvar='Omegatc';  % nan
% inoutvar='qsc';  % 1e-12
% inoutvar='qst';  % 1e-12
% inoutvar='term3';
% inoutvar='qs';


% select qtrans model to be tested
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% Use one of the duck94 test case at time step 'targetn' as background
% state.  Note, model input data are pre-cached as mat-file.
duck94Case='c'; targetn=250;
% duck94Case='b'; targetn=200;
load(['../../../../test_cases/duck94/obsdataCache/obsdata_case' duck94Case '.mat']);
modelinput=initModelInputs(duck94Case,grid,sedmodel);

% extract relevant model inputs for time step 'targetn'
dt=diff([hydroobs(targetn+[0 1]).dnum_est])*24*3600;
nsubsteps=1;
x=grid.x;
nx=length(x);
H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,hydroobs(targetn).dnum_est);
theta0=interp1(waves8m.dnum_est,waves8m.theta0,hydroobs(targetn).dnum_est);
omega =interp1(waves8m.dnum_est,waves8m.sigmam,hydroobs(targetn).dnum_est);
tau_wind=interp1(windEOP.dnum_est,windEOP.tau,hydroobs(targetn).dnum_est);
tau_wind=repmat(tau_wind,nx,1);
dgamma=ones(nx,1)*.01;
dAw=ones(nx,1)*.01;
dSw=ones(nx,1)*.01;
detady=ones(nx,1)*.001;
tide=interp1(waves8m.dnum_est,waves8m.tide,hydroobs(targetn).dnum_est);
d50      =modelinput.d50      ;
d90      =modelinput.d90      ;
ka_drag  =modelinput.ka_drag  ;
beta0    =modelinput.beta0    ;
betaType =modelinput.betaType ;
gammaType=modelinput.gammaType;
param=modelinput.params;
h=grid.h+tide;
h(h<.5)=.5;  % min depth constraint

% OPTIONAL: override default hydro model formulations to test their TL-AD
% gammaType=2003;
% betaType='const';
% param.streamingType='v';  % choose either 'n' or 'v'

% run hydro models to get inputs
[Hrms,theta,vbar,kabs,Ew,Er,Dr,hydro_bkgd]=hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,beta0,gammaType,betaType);
Hmo=1.4*Hrms;
[Aw,Sw,Uw,uwave_bkgd]=Uwave_ruessink2012_params(Hmo,kabs,omega,h);
rho=1030;
c=omega./kabs;
kvec=cat(2,kabs(:).*cos(theta(:)),kabs(:).*sin(theta(:)));
ubarx=-(Ew+2*Er)./(rho*c.*h);  % e.g., Dubarbier et al. (2015) eqn 8
ubar(:,1)=ubarx;
ubar(:,2)=vbar;
delta=ones(nx,1)*.2;  % init
udelta=zeros(nx,2);   % init
for i=1:nx
  if(Dr(i)>0)
    [udelta(i,:),delta(i),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),kvec(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           tau_wind(i,:),Dr(i),param.fv,param.ks,d50);
  end
end
tanbeta=calcTanbeta(x,h)';
ws=ws_brownLawler(.8*d50);

% bkgd NL model run
[Q,bkgd]=qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(13*nx+7,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  tl_d50        =F(0*nx+0+[1:nx],i);
  tl_d90        =F(1*nx+0+[1:nx],i);
  tl_h          =F(2*nx+0+[1:nx],i);
  tl_tanbeta    =F(3*nx+0+[1:nx],i);
  tl_Hrms       =F(4*nx+0+[1:nx],i);
  tl_kabs       =F(5*nx+0+[1:nx],i);
  tl_omega      =F(6*nx+1       ,i);
  tl_udelta(:,1)=F(6*nx+1+[1:nx],i);
  tl_udelta(:,2)=F(7*nx+1+[1:nx],i);
  tl_delta      =F(8*nx+1+[1:nx],i);
  tl_ws         =F(9*nx+1+[1:nx],i);
  tl_Aw         =F(10*nx+1+[1:nx],i);
  tl_Sw         =F(11*nx+1+[1:nx],i);
  tl_Uw         =F(12*nx+1+[1:nx],i);
  tl_param.n    =F(13*nx+2       ,i);
  tl_param.m    =F(13*nx+3       ,i);
  tl_param.xi   =F(13*nx+4       ,i);
  tl_param.alpha=F(13*nx+5       ,i);
  if(isfield(param,'Cc'))
    tl_param.Cc   =F(13*nx+6       ,i);
    tl_param.Cf   =F(13*nx+7       ,i);
  end

  % TL model: TL*F
  tl_qs=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_delta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_d50,ad_d90,ad_h,ad_tanbeta,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_delta,ad_ws,ad_Aw,ad_Sw,ad_Uw,ad_param] = ...
      ad_qtrans_vanderA(tl_qs,bkgd);%,inoutvar);

  % create output vector
  g(0*nx+0+[1:nx],i) =ad_d50        ;
  g(1*nx+0+[1:nx],i) =ad_d90        ;
  g(2*nx+0+[1:nx],i) =ad_h          ;
  g(3*nx+0+[1:nx],i) =ad_tanbeta    ;
  g(4*nx+0+[1:nx],i) =ad_Hrms       ;
  g(5*nx+0+[1:nx],i) =ad_kabs       ;
  g(6*nx+1       ,i) =ad_omega      ;
  g(6*nx+1+[1:nx],i) =ad_udelta(:,1);
  g(7*nx+1+[1:nx],i) =ad_udelta(:,2);
  g(8*nx+1+[1:nx],i) =ad_delta      ;
  g(9*nx+1+[1:nx],i) =ad_ws         ;
  g(10*nx+1+[1:nx],i) =ad_Aw         ;
  g(11*nx+1+[1:nx],i) =ad_Sw         ;
  g(12*nx+1+[1:nx],i)=ad_Uw         ;
  g(13*nx+2       ,i)=ad_param.n    ;
  g(13*nx+3       ,i)=ad_param.m    ;
  g(13*nx+4       ,i)=ad_param.xi   ;
  g(13*nx+5       ,i)=ad_param.alpha;
  if(isfield(param,'Cc'))
    g(13*nx+6       ,i)=ad_param.Cc   ;
    g(13*nx+7       ,i)=ad_param.Cf   ;
  end

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
