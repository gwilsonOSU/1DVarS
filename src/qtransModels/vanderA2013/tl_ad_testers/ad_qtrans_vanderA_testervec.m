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
for i=1:nx
  if(Dr(i)==0)
    udelta(i,:)=[0 0];
  else
    [udelta(i,:),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),kvec(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           windW(i,:),Dr(i),param.fv,param.ks,d50);
  end
end

% bkgd NL model run
param.n=1.2;
param.m=11;
param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
param.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)
[Q,bkgd]=qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta,ws,Aw,Sw,Uw,param);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(11*nx+5,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  tl_d50        =F(0*nx+0+[1:nx],i);
  tl_d90        =F(1*nx+0+[1:nx],i);
  tl_h          =F(2*nx+0+[1:nx],i);
  tl_Hrms       =F(3*nx+0+[1:nx],i);
  tl_kabs       =F(4*nx+0+[1:nx],i);
  tl_omega      =F(5*nx+1       ,i);
  tl_udelta(:,1)=F(5*nx+1+[1:nx],i);
  tl_udelta(:,2)=F(6*nx+1+[1:nx],i);
  tl_ws         =F(7*nx+1+[1:nx],i);
  tl_Aw         =F(8*nx+1+[1:nx],i);
  tl_Sw         =F(9*nx+1+[1:nx],i);
  tl_Uw         =F(10*nx+1+[1:nx],i);
  tl_param.n    =F(11*nx+2       ,i);
  tl_param.m    =F(11*nx+3       ,i);
  tl_param.xi   =F(11*nx+4       ,i);
  tl_param.alpha=F(11*nx+5       ,i);

  % TL model: TL*F
  tl_qs=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_d50,ad_d90,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_Aw,ad_Sw,ad_Uw,ad_param] = ...
      ad_qtrans_vanderA(tl_qs,bkgd);%,inoutvar);

  % create output vector
  g(0*nx+0+[1:nx],i) =ad_d50        ;
  g(1*nx+0+[1:nx],i) =ad_d90        ;
  g(2*nx+0+[1:nx],i) =ad_h          ;
  g(3*nx+0+[1:nx],i) =ad_Hrms       ;
  g(4*nx+0+[1:nx],i) =ad_kabs       ;
  g(5*nx+1       ,i) =ad_omega      ;
  g(5*nx+1+[1:nx],i) =ad_udelta(:,1);
  g(6*nx+1+[1:nx],i) =ad_udelta(:,2);
  g(7*nx+1+[1:nx],i) =ad_ws         ;
  g(8*nx+1+[1:nx],i) =ad_Aw         ;
  g(9*nx+1+[1:nx],i) =ad_Sw         ;
  g(10*nx+1+[1:nx],i)=ad_Uw         ;
  g(11*nx+2       ,i)=ad_param.n    ;
  g(11*nx+3       ,i)=ad_param.m    ;
  g(11*nx+4       ,i)=ad_param.xi   ;
  g(11*nx+5       ,i)=ad_param.alpha;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
