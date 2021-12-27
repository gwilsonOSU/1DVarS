%
% symmetry test
%
clear
addpath(genpath('../../..'))

% % TEST: select an i/o variable for testing
% inoutvar='Hmo';
% inoutvar='Uw';
% inoutvar='utilde';
% inoutvar='uH';
% inoutvar='Au_nums';
% inoutvar='Au_dens_1';
% inoutvar='Au_dens';
% inoutvar='Au';
% inoutvar='Aw';
% inoutvar='utot';
% inoutvar='qa';
% % inoutvar='utabs';
% % inoutvar='utotabs';
% % inoutvar='qb1';
% % inoutvar='qb2';
% % inoutvar='qb3';
% inoutvar='qb';
% % inoutvar='qs1';
% % inoutvar='qs2';
% % inoutvar='qs3';
% inoutvar='qs';
% inoutvar='q';
% inoutvar='Q';

load ~/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

param.Cw=0.00483 ;
param.Cc=0.02002 ;
param.Cf=0.01173 ;
param.Ka=0.631e-4;
param.dt=60*60;  % one hour time step

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

% background NL model run
[Q,Qb,Qs,Qa,bkgd] = ...
    qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta,ws,Aw,Sw,Uw,...
                     param.Cw,param.Cc,param.Cf,param.Ka);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(10*nx+6,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_tanbeta=F(0*nx+[1:nx],i);
  tl_h      =F(1*nx+[1:nx],i);
  tl_Hrms   =F(2*nx+[1:nx],i);
  tl_kabs   =F(3*nx+[1:nx],i);
  tl_omega  =F(4*nx+1,i);
  tl_udelta =reshape(F(4*nx+1+[1:(2*nx)],i),[nx 2]);
  tl_ws     =F(6*nx+2+[1:nx],i);
  tl_Aw     =F(7*nx+2+[1:nx],i);
  tl_Sw     =F(8*nx+2+[1:nx],i);
  tl_Uw     =F(9*nx+2+[1:nx],i);
  tl_param.Cw=F(10*nx+3,i);
  tl_param.Cc=F(10*nx+4,i);
  tl_param.Cf=F(10*nx+5,i);
  tl_param.Ka=F(10*nx+6,i);
  tl_Q =tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,...
                            tl_param.Cw,tl_param.Cc,tl_param.Cf,tl_param.Ka,...
                            bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_Aw,ad_Sw,ad_Uw,...
   ad_Cw,ad_Cc,ad_Cf,ad_Ka] = ...
      ad_qtrans_dubarbier(tl_Q,bkgd);%,inoutvar);
  % create output vector
  g(0*nx+[1:nx],i)    =ad_tanbeta;
  g(1*nx+[1:nx],i)    =ad_h      ;
  g(2*nx+[1:nx],i)    =ad_Hrms   ;
  g(3*nx+[1:nx],i)    =ad_kabs   ;
  g(4*nx+1,i)         =ad_omega  ;
  g(4*nx+1+[1:(2*nx)],i)=ad_udelta(:);
  g(6*nx+2+[1:nx],i)  =ad_ws     ;
  g(7*nx+2+[1:nx],i)  =ad_Aw     ;
  g(8*nx+2+[1:nx],i)  =ad_Sw     ;
  g(9*nx+2+[1:nx],i)  =ad_Uw     ;
  g(10*nx+3,i)         =ad_Cw;
  g(10*nx+4,i)         =ad_Cc;
  g(10*nx+5,i)         =ad_Cf;
  g(10*nx+6,i)         =ad_Ka;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
