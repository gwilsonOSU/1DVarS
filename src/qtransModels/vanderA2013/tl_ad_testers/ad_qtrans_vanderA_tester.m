%
% TL-AD symmetry test
%
clear

% % TEST-CODE: choose a test variable.  These are in the order they appear in
% % TL model.  This is a very effective debugging technique for the AD model.
% % Requires having an auxiliary input for the TL and AD codes that override
% % the default output variable.  See {tl,ad}_qtrans_vanderA.m 'TEST-CODE'
% % blocks (commented out for production code)
% inoutvar='Hmo';
% inoutvar='c';
% inoutvar='uw2mean';
% inoutvar='uhat';
% inoutvar='asinarg ';
% inoutvar='Tc';
% inoutvar='Tt';
% inoutvar='T';  % nan
% inoutvar='Tcu';
% inoutvar='Ttu';
% inoutvar='uhatc';
% inoutvar='uhatt';
% inoutvar='Dstar';
% inoutvar='theta_cr';
% inoutvar='ahat';
% inoutvar='utildecr';
% inoutvar='utildetr';
% inoutvar='psihatc';
% inoutvar='psihatt';
% inoutvar='nlambda';  % nan
% inoutvar='eta';  % nan
% inoutvar='lambda';  % nan
% inoutvar='udabs';
% inoutvar='alpha';
% inoutvar='fwdc';
% inoutvar='fwdt';
% inoutvar='ucrvec';
% inoutvar='utrvec';
% inoutvar='ucrabs';
% inoutvar='utrabs';
% inoutvar='thetac';
% inoutvar='thetat';
% inoutvar='thetahatc';
% inoutvar='thetahatt';
% inoutvar='Pc';
% inoutvar='Pt';
% inoutvar='absthetac';
% inoutvar='absthetat';
% inoutvar='term1';
% inoutvar='term2';
% inoutvar='term3';
% inoutvar='qs';

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
h=waves.h(i);
Hrms=waves.H(i);
ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
ubar(2)=waves.v(i);
kabs=waves.k(i);
omega=waves.sigma;
detady=.001;
windW =[1 1]*.5;
Dr=waves.eps_r(i)*1030;
fv=1;
d50=180e-6;
d90=400e-6;
k=kabs*[cos(waves.theta(i)) sin(waves.theta(i))];
udelta=udelta_reniers2004(ubar,k,omega,h,...
                               Hrms,detady,windW,...
                               Dr,fv,d50);
ws=ws_brownLawler(d50);
param.n=1.2;
param.m=11;
param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
param.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)
[q,bkgd]=qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta,ws,param);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F = eps*rand(13,n);  % 1st dim is number of tl input parameters
for i=1:n
  % TL model: TL*F
  tl_d50        =F(1,i);
  tl_d90        =F(2,i);
  tl_h          =F(3,i);
  tl_Hrms       =F(4,i);
  tl_kabs       =F(5,i);
  tl_omega      =F(6,i);
  tl_udelta     =F(7:8,i)';
  tl_ws         =F(9,i);
  tl_param.n    =F(10,i);
  tl_param.m    =F(11,i);
  tl_param.xi   =F(12,i);
  tl_param.alpha=F(13,i);
  tl_qs=tl_qtrans_vanderA(tl_d50,...
                                    tl_d90,...
                                    tl_h,...
                                    tl_Hrms,...
                                    tl_kabs,...
                                    tl_omega,...
                                    tl_udelta,...
                                    tl_ws,...
                                    tl_param,...
                                    bkgd);%,inoutvar);
  % AD model: g=AD*(TL*F)
  [ad_d50,ad_d90,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_param] = ...
      ad_qtrans_vanderA(tl_qs,bkgd);%,inoutvar);
  g(:,i)=[ad_d50;ad_d90;ad_h;ad_Hrms;ad_kabs;ad_omega;ad_udelta(:);ad_ws;
          ad_param.n;ad_param.m;ad_param.xi;ad_param.alpha];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
