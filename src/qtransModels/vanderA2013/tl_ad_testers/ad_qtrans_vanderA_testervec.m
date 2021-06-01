%
% symmetry test
%
clear

% % TEST-CODE: choose a test variable.  These are in the order they appear in
% % TL model.  This is a very effective debugging technique for the AD model.
% % Requires having an auxiliary input for the TL and AD codes that override
% % the default output variable.  See {tl,ad}_qtrans_vanderA.m 'TEST-CODE'
% % blocks (commented out for production code)
% inoutvar='d50';
% inoutvar='d90';
% inoutvar='h';
% inoutvar='Hrms';
% inoutvar='kabs';
% inoutvar='omega';
% inoutvar='udelta';
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
% inoutvar='tauwRe';
% inoutvar='streamingEffect';
% inoutvar='thetacx';
% inoutvar='thetatx';
% inoutvar='thetahatc';
% inoutvar='thetahatt';
% inoutvar='deltasc';
% inoutvar='deltast';
% inoutvar='deltasc';
% inoutvar='deltast';
% inoutvar='b';
% inoutvar='RR';
% inoutvar='worb1c';
% inoutvar='worb1t';
% inoutvar='worb2c';
% inoutvar='worb2t';
% inoutvar='t1ca';
% inoutvar='t1c';
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
% inoutvar='qsc';
% inoutvar='qst';
% inoutvar='term3';
% inoutvar='qs';

load ~/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% define bkgd variables
i=50:75;  % gridpoints
d50=180e-6;
d90=400e-6;
x=waves.x(i);
h=filtfilt(ones(5,1)/5,1,waves.h(i));
tanbeta=calcTanbeta(x,h);
Hrms=waves.H(i);
Ew=waves.E(i)*rho;
Er=waves.Er(i)*rho;
Dr=waves.eps_r(i)*rho;
ubar(:,1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
ubar(:,2)=waves.v(i);
omega=waves.sigma;
kabs=waves.k(i);
kvec(:,1)=waves.k(i).*cos(waves.theta(i));
kvec(:,2)=waves.k(i).*sin(waves.theta(i));
theta=waves.theta(i);
windW=zeros(length(x),2);
detady=zeros(length(x),2);
ws=ws_brownLawler(d50);

% reniers model for udelta
nx=length(x);
param.fv=.1;
for i=1:nx
  if(Dr(i)==0)
    udelta(i,:)=[0 0];
  else
    [udelta(i,:),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),kvec(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           windW(i,:),Dr(i),param.fv,d50);
  end
end

% bkgd NL model run
param.n=1.2;
param.m=11;
param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
param.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)
for j=1:nx
  [Q(j),bkgd(j)]=qtrans_vanderA(d50,d90,h(j),Hrms(j),kabs(j),omega,udelta(j,:),ws,param);
end

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(5*nx+8,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_h      =F(0*nx+[1:nx],i);
  tl_Hrms   =F(1*nx+[1:nx],i);
  tl_kabs   =F(2*nx+[1:nx],i);
  tl_udelta =reshape(F(3*nx+[1:(2*nx)],i),[nx 2]);
  tl_ws     =F(5*nx+1,i);
  tl_param.n=F(5*nx+2,i);
  tl_param.m=F(5*nx+3,i);
  tl_param.xi=F(5*nx+4,i);
  tl_param.alpha=F(5*nx+5,i);
  tl_d50=F(5*nx+6,i);
  tl_d90=F(5*nx+7,i);
  tl_omega=F(5*nx+8,i);
  for j=1:nx
    tl_Q(j)=tl_qtrans_vanderA(tl_d50,tl_d90,tl_h(j),tl_Hrms(j),...
                                           tl_kabs(j),tl_omega,tl_udelta(j,:),tl_ws,tl_param,...
                                           bkgd(j));%,inoutvar);
  end

  % AD model: g=AD*(TL*F)
  ad_d50=0;
  ad_d90=0;
  ad_h=zeros(nx,1);
  ad_Hrms=zeros(nx,1);
  ad_kabs=zeros(nx,1);
  ad_udelta=zeros(nx,2);
  ad_ws=0;
  ad_param.n    =0;
  ad_param.m    =0;
  ad_param.xi   =0;
  ad_param.alpha=0;
  ad_omega=0;
  for j=nx:-1:1
    [ad1_d50,ad1_d90,ad1_h,ad1_Hrms,ad1_kabs,ad1_omega,ad1_udelta,ad1_ws,ad1_param] = ...
        ad_qtrans_vanderA(tl_Q(j),bkgd(j));%,inoutvar);
    ad_d50        =ad_d50        +ad1_d50        ;
    ad_d90        =ad_d90        +ad1_d90        ;
    ad_h(j)       =ad_h(j)       +ad1_h          ;
    ad_Hrms(j)    =ad_Hrms(j)    +ad1_Hrms       ;
    ad_kabs(j)    =ad_kabs(j)    +ad1_kabs       ;
    ad_udelta(j,:)=ad_udelta(j,:)+ad1_udelta     ;
    ad_ws         =ad_ws         +ad1_ws         ;
    ad_param.n    =ad_param.n    +ad1_param.n    ;
    ad_param.m    =ad_param.m    +ad1_param.m    ;
    ad_param.xi   =ad_param.xi   +ad1_param.xi   ;
    ad_param.alpha=ad_param.alpha+ad1_param.alpha;
    ad_omega=ad_omega+ad1_omega;
  end
  % create output vector
  g(0*nx+[1:nx],i)    =ad_h          ;
  g(1*nx+[1:nx],i)    =ad_Hrms       ;
  g(2*nx+[1:nx],i)    =ad_kabs       ;
  g(3*nx+[1:(2*nx)],i)=ad_udelta(:)  ;
  g(5*nx+1,i)         =ad_ws         ;
  g(5*nx+2,i)         =ad_param.n    ;
  g(5*nx+3,i)         =ad_param.m    ;
  g(5*nx+4,i)         =ad_param.xi   ;
  g(5*nx+5,i)         =ad_param.alpha;
  g(5*nx+6,i)         =ad_d50;
  g(5*nx+7,i)         =ad_d90;
  g(5*nx+8,i)         =ad_omega;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
