%
% symmetry test
%
clear

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

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

params.Cw=0.00483 ;
params.Cc=0.02002 ;
params.Cf=0.01173 ;
params.Ka=0.631e-4;
params.dt=60*60;  % one hour time step
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)

% define bkgd variables
i=1:20;  % gridpoints
d50=250e-6;
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
for i=1:nx
  if(Dr(i)==0)
    udelta(i,:)=[0 0];
  else
    [udelta(i,:),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),kvec(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           windW(i,:),Dr(i),params.fv,d50);
  end
end

% background NL model run
[Q,Qb,Qs,Qa,bkgd] = ...
    qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta,ws,...
                     params.Cw,params.Cc,params.Cf,params.Ka);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(6*nx+6,n);  % 1st dim is number of tl input parameters
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
  tl_ws     =F(6*nx+2,i);
  tl_params.Cw=F(6*nx+3,i);
  tl_params.Cc=F(6*nx+4,i);
  tl_params.Cf=F(6*nx+5,i);
  tl_params.Ka=F(6*nx+6,i);
  tl_Q =tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,...
                            tl_params.Cw,tl_params.Cc,tl_params.Cf,tl_params.Ka,...
                            bkgd);%,inoutvar);

  % AD model: g=AD*(TL*F)
  [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,...
   ad_Cw,ad_Cc,ad_Cf,ad_Ka] = ...
      ad_qtrans_dubarbier(tl_Q,bkgd);%,inoutvar);
  % create output vector
  g(0*nx+[1:nx],i)    =ad_tanbeta;
  g(1*nx+[1:nx],i)    =ad_h      ;
  g(2*nx+[1:nx],i)    =ad_Hrms   ;
  g(3*nx+[1:nx],i)    =ad_kabs   ;
  g(4*nx+1,i)         =ad_omega  ;
  g(4*nx+1+[1:(2*nx)],i)=ad_udelta(:);
  g(6*nx+2,i)         =ad_ws     ;
  g(6*nx+3,i)         =ad_Cw;
  g(6*nx+4,i)         =ad_Cc;
  g(6*nx+5,i)         =ad_Cf;
  g(6*nx+6,i)         =ad_Ka;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
