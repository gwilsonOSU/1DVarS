%
% perturbation test of tl_qtrans_vanderA.m
%
clear
addpath(genpath('../../..'));

% TEST: choose output variable to test
outvarlist={};
% outvarlist{end+1}='phi_r2012';
% outvarlist{end+1}='r_r2012';
% % outvarlist{end+1}='uhat';
% % outvarlist{end+1}='T';
% outvarlist{end+1}='phidc';
% outvarlist{end+1}='phiuc';
% outvarlist{end+1}='tdc';  % sometimes off
% outvarlist{end+1}='tuc';
% outvarlist{end+1}='Tc';
% outvarlist{end+1}='Tt';
% outvarlist{end+1}='Tcu';
% outvarlist{end+1}='Ttu';
% outvarlist{end+1}='uhatc';
% outvarlist{end+1}='uhatt';
% outvarlist{end+1}='Dstar';
% outvarlist{end+1}='theta_cr';
% outvarlist{end+1}='ahat';  % small bias
% outvarlist{end+1}='utildecr';
% outvarlist{end+1}='utildetr';
% outvarlist{end+1}='psihatc';  % small bias
% outvarlist{end+1}='psihatt';  % small bias
% outvarlist{end+1}='nlambda';
% outvarlist{end+1}='eta';
% outvarlist{end+1}='lambda';
% outvarlist{end+1}='alpha';
% outvarlist{end+1}='ksd';
% outvarlist{end+1}='fwc';
% outvarlist{end+1}='fwdc';
% outvarlist{end+1}='fwdt';
% % outvarlist{end+1}='ucrvec';
% % outvarlist{end+1}='utrvec';
% outvarlist{end+1}='ucrabs';
% outvarlist{end+1}='utrabs';
% outvarlist{end+1}='thetac';
% outvarlist{end+1}='thetat';
% outvarlist{end+1}='streamingEffect';
% outvarlist{end+1}='thetacx';
% outvarlist{end+1}='thetatx';
% outvarlist{end+1}='thetahatc';
% outvarlist{end+1}='thetahatt';
% outvarlist{end+1}='b';
% outvarlist{end+1}='RR';
% outvarlist{end+1}='deltasc';
% outvarlist{end+1}='etawc';
% outvarlist{end+1}='worb1c';
% outvarlist{end+1}='worb1t';
% outvarlist{end+1}='worb2c';
% outvarlist{end+1}='worb2t';
% outvarlist{end+1}='t1ca';
% outvarlist{end+1}='t1ta';
% outvarlist{end+1}='t1cb';
% outvarlist{end+1}='t1tb';
% outvarlist{end+1}='t1c';
% outvarlist{end+1}='t1t';
% outvarlist{end+1}='t2cb';
% outvarlist{end+1}='t2tb';
% outvarlist{end+1}='t2ca';
% outvarlist{end+1}='t2ta';
% outvarlist{end+1}='t2c';
% outvarlist{end+1}='t2t';
% outvarlist{end+1}='worbc';
% outvarlist{end+1}='worbt';
% outvarlist{end+1}='wsc';
% outvarlist{end+1}='wst';
% outvarlist{end+1}='Pc';
% outvarlist{end+1}='Pt';
% outvarlist{end+1}='qsc';
% outvarlist{end+1}='qst';
% outvarlist{end+1}='qsCc';
% outvarlist{end+1}='qsCf';
% outvarlist{end+1}='qs';
% for ii=1:length(outvarlist)
%   outvar=outvarlist{ii};
%   for jj=1
%   disp(outvar)

% select qtrans model to be tested
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% Use one of the duck94 test case at time step 'targetn' as background
% state.  Note, model input data are pre-cached as mat-file.
% duck94Case='c'; targetn=250;
duck94Case='b'; targetn=200;
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
dAw=ones(nx,1)*01;
dSw=ones(nx,1)*01;
detady=ones(nx,1)*.0001;
tide=interp1(waves8m.dnum_est,waves8m.tide,hydroobs(targetn).dnum_est);
d50      =modelinput.d50      ;
d90      =modelinput.d90      ;
ka_drag  =modelinput.ka_drag  ;
beta0    =modelinput.beta0    ;
betaType =modelinput.betaType ;
gammaType=modelinput.gammaType;
param=modelinput.params;
h=grid.h+tide;
h(h<.75)=.75;  % min depth constraint

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
[q,bkgd]=qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);%,outvar);

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_h          =[bkgd.h      ]'    *frac_tl*myrand();
tl_d50        =[bkgd.d50    ]'    *frac_tl*myrand();
tl_d90        =[bkgd.d90    ]'    *frac_tl*myrand();
tl_h          =[bkgd.h      ]'    *frac_tl*myrand();
tl_tanbeta    =[bkgd.tanbeta]'    *frac_tl*myrand();
tl_Hrms       =[bkgd.Hrms   ]'    *frac_tl*myrand();
tl_kabs       =[bkgd.kabs   ]'    *frac_tl*myrand();
tl_omega      =bkgd(1).omega      *frac_tl*myrand();
tl_udelta     =reshape([bkgd.udelta ]'    *frac_tl*myrand(),[2 nx])';
tl_delta      =[bkgd.delta  ]'    *frac_tl*myrand();
tl_ws         =[bkgd.ws     ]'    *frac_tl*myrand();
tl_Aw         =[bkgd.Aw     ]'    *frac_tl*myrand();
tl_Sw         =[bkgd.Sw     ]'    *frac_tl*myrand();
tl_Uw         =[bkgd.Uw     ]'    *frac_tl*myrand();
pp=bkgd(1).param;
tl_param.fv   =[pp.fv   ]*frac_tl*myrand();
tl_param.n    =[pp.n    ]*frac_tl*myrand();
tl_param.m    =[pp.m    ]*frac_tl*myrand();
tl_param.xi   =[pp.xi   ]*frac_tl*myrand();
tl_param.alpha=[pp.alpha]*frac_tl*myrand();
if(isfield(param,'Cc'))
  tl_param.Cc   =[pp.Cc   ]*frac_tl*myrand();
  tl_param.Cf   =[pp.Cf   ]*frac_tl*myrand();
end

% pack perturbed param into struct for NL model
param_prime.fv   =param.fv   +tl_param.fv   ;
param_prime.n    =param.n    +tl_param.n    ;
param_prime.m    =param.m    +tl_param.m    ;
param_prime.xi   =param.xi   +tl_param.xi   ;
param_prime.alpha=param.alpha+tl_param.alpha;
if(isfield(param,'Cc'))
  param_prime.Cc   =param.Cc   +tl_param.Cc   ;
  param_prime.Cf   =param.Cf   +tl_param.Cf   ;
end

% run NL model with perturbation
[q_prime,bkgd_prime] = ...
    qtrans_vanderA(d50+tl_d50,d90+tl_d90,h+tl_h,tanbeta+tl_tanbeta,Hrms+tl_Hrms,...
                   kabs+tl_kabs,omega+tl_omega,udelta+tl_udelta,delta+tl_delta,...
                   ws+tl_ws,Aw+tl_Aw,Sw+tl_Sw,Uw+tl_Uw,param_prime);%,outvar);
tl_q_true = q_prime - q;

% run TL model
tl_q = tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,...
                         tl_omega,tl_udelta,tl_delta,...
                         tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd);%,outvar);

% compare scatter
imax=150;  % ignore shoreline
clf
subplot(121)
plot(tl_q_true(1:imax),tl_q(1:imax),'.')
hold on
axis equal tight
ax=axis;
plot(ax([1 2]),ax([1 2]),'k--')
hold off
xlabel('true')
ylabel('predicted')

% compare vector
subplot(122), hold on
plot(tl_q_true)
plot(tl_q)
legend('true','predicted')
xlim([1 imax])

% drawnow;
% pause(.1);
% end
% end
