%
% perturbation test of tl_qtrans_vanderA.m
%
clear
addpath(genpath('../../..'));

% TEST: choose output variable to test
% outvar='Hmo';
% outvar='c';
% outvar='Uw';
% outvar='phi_r2012';
% outvar='r_r2012';  % off, but close.  Analytical derivs seem off
% % outvar='uhat';
% % outvar='T';
% outvar='Tc';  %fixed, ok
% outvar='Tt';  %Ok
% outvar='Tcu';  %ok
% outvar='Ttu'; %ok
% % outvar='uhatc'; % depends on r_r2012
% % outvar='uhatt'; % depends on r_r2012
% % outvar='Dstar';
% %outvar='theta_cr'; % depends on r_r2012
% outvar='ahat'; %ok
% % outvar='utildecr'; % depends on r_r2012
% % outvar='utildetr'; % depends on r_r2012
% outvar='psihatc'; %ok
% outvar='psihatt'; %ok
% % outvar='nlambda';
% % outvar='eta';
% % outvar='lambda';
% % outvar='alpha';
% outvar='fwdc';  % ok but regime dependent
% outvar='fwdt';  % ok but regime dependent
% % outvar='ucrvec';  % can't test nonscalars
% % outvar='utrvec';  % can't test nonscalars
% outvar='ucrabs'; %ok
% outvar='utrabs'; %ok
% outvar='thetac';  % ok but regime dependent
% outvar='thetat';  % ok but regime dependent
% outvar='streamingEffect';
% outvar='thetacx';
% outvar='thetatx';
% outvar='thetahatc';
% % outvar='thetahatt';
% % outvar='b';
% % outvar='RR';
% % outvar='worb1c';
% % outvar='worb1t';
% % outvar='worb2c';
% % outvar='worb2t';
% % outvar='t1ca';
% % outvar='t1ta';
% % outvar='t1cb';
% % outvar='t1tb';
% % outvar='t1c';
% % outvar='t1t';
% % outvar='t2cb';
% % outvar='t2tb';
% % outvar='t2ca';
% % outvar='t2ta';
% % outvar='t2c';
% % outvar='t2t';
% % outvar='worbc';
% % outvar='worbt';
% % outvar='wsc';
% % outvar='wst';
% % outvar='Pc';
% % outvar='Pt';
% % outvar='absthetac';
% % outvar='absthetat';
% % outvar='qsc';  % WAY off
% % outvar='qst';  % way off
% % outvar='term3';
% % outvar='uwmo';
% % outvar='qs2';  %ok
% % outvar='qs3';  %ok
% % outvar='qsCc'; %ok
% % outvar='qsCf'; %Ok
% outvar='qs';

param.streamingType='v';  % choose either 'n' or 'v'.  NOTE, should test both

load ~/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define bkgd variables
x=waves.x;
nx=length(x);
d50=180e-6*ones(nx,1);
d90=400e-6*ones(nx,1);
h=filtfilt(ones(5,1)/5,1,waves.h);
tanbeta=calcTanbeta(x,h);
Hrms=waves.H;
Ew=waves.E*rho;
Er=waves.Er*rho;
Dr=waves.eps_r*rho;
ubar(:,1)=-(waves.E+2*waves.Er)./(waves.c.*waves.h);
ubar(:,2)=waves.v;
omega=waves.sigma;
kabs=waves.k;
kvec(:,1)=waves.k.*cos(waves.theta);
kvec(:,2)=waves.k.*sin(waves.theta);
theta=waves.theta;
windW=zeros(length(x),2);
detady=zeros(length(x),2);
ws=ws_brownLawler(d50);
kabs=sqrt(sum(kvec.^2,2));
[Aw,Sw,Uw,uwave_wksp]=Uwave_ruessink2012_params(1.4*Hrms,kabs,omega,h);

% reniers model for udelta
nx=length(x);
param.fv=.1;
param.ks=.0083;
delta=ones(nx,1)*.2;  % init
udelta=zeros(nx,2);   % init
for i=1:nx
  if(Dr(i)>0)
    [udelta(i,:),delta(i),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),kvec(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           windW(i,:),Dr(i),param.fv,param.ks,d50);
  end
end
tanbeta=calcTanbeta(x,h)';

% bkgd NL model run
param.n=1.2;
param.m=11;
param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
param.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)
% param.Cc=0.01;  % stirring+undertow effect.  Comment out to use auto-above-WBL transport
% param.Cf=0.01;  % stirring+slope effect
[q,bkgd]=qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);%,outvar);

% choose reasonable perturbations
frac_tl = 0.001;
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

% pack perturbed params into struct for NL model
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
clf
subplot(121)
plot(tl_q_true,tl_q,'.')
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
