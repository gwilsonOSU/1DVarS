%
% TL-AD symmetry test
%
clear

load ~/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
bkgd.h=waves.h(i);
bkgd.Hrms=waves.H(i);
bkgd.ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
bkgd.ubar(2)=waves.v(i);
bkgd.omega=waves.sigma;
bkgd.detady=.001;
bkgd.windW =[1 1]*.5;
bkgd.Dr=waves.eps_r(i)*1030;
bkgd.fv=1;
bkgd.d50=180e-6;
bkgd.d90=400e-6;
bkgd.kabs=waves.k(i);
bkgd.k=waves.k(i)*[cos(waves.theta(i)) sin(waves.theta(i))];
bkgd.udelta=udelta_reniers2004(bkgd.ubar,bkgd.k,bkgd.omega,bkgd.h,...
                               bkgd.Hrms,bkgd.detady,bkgd.windW,...
                               bkgd.Dr,bkgd.fv,bkgd.d50);
bkgd.ws=ws_brownLawler(bkgd.d50);
bkgd.param.n=1.2;
bkgd.param.m=11;
bkgd.param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
bkgd.param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19

% bkgd state.  Need to run the "full" model, can't go straight to the
% shields solver submodel (vanderA_shields.m)
[qs,bkgd_qt]=qtrans_vanderA(bkgd.d50,...
                            bkgd.d90,...
                            bkgd.h,...
                            bkgd.Hrms,...
                            bkgd.kabs,...
                            bkgd.omega,...
                            bkgd.udelta,...
                            bkgd.ws,...
                            bkgd.param);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F = eps*rand(8,n);  % 1st dim is number of tl input parameters
for i=1:n
  % perturbations
  tl_d50   =F(1,i);
  tl_d90   =F(2,i);
  tl_udelta=F(3,i);
  tl_uhat  =F(4,i);
  tl_mu    =F(5,i);
  tl_eta   =F(6,i);
  tl_lambda=F(7,i);
  tl_ahat  =F(8,i);
  % TL model: TL*F
  [tl_theta_av,tl_ksd,tl_ksw,tl_fd,tl_fw] = ...
      tl_vanderA_shields(tl_d50,tl_d90,tl_udelta,tl_uhat,...
                         tl_mu,tl_eta,tl_lambda,tl_ahat,...
                         bkgd_qt.d50,...
                         bkgd_qt.d90,...
                         bkgd_qt.udabs,...
                         bkgd_qt.uhat,...
                         bkgd_qt.delta,...
                         bkgd_qt.mu,...
                         bkgd_qt.eta,...
                         bkgd_qt.lambda,...
                         bkgd_qt.ahat,...
                         bkgd_qt.ksd,...
                         bkgd_qt.ksw,...
                         bkgd_qt.fd,...
                         bkgd_qt.fw,...
                         bkgd_qt.branch_A1,...
                         bkgd_qt.branch_A4,...
                         bkgd_qt.branch_A5);
  % AD model: g=AD*(TL*F)
  [ad_d50,ad_d90,ad_udelta,ad_uhat,...
   ad_mu,ad_eta,ad_lambda,ad_ahat] = ...
      ad_vanderA_shields(tl_theta_av,tl_ksd,tl_ksw,tl_fd,tl_fw,...
                         bkgd_qt.d50,...
                         bkgd_qt.udabs,...
                         bkgd_qt.uhat,...
                         bkgd_qt.delta,...
                         bkgd_qt.mu,...
                         bkgd_qt.eta,...
                         bkgd_qt.lambda,...
                         bkgd_qt.ahat,...
                         bkgd_qt.ksd,...
                         bkgd_qt.ksw,...
                         bkgd_qt.fd,...
                         bkgd_qt.fw,...
                         bkgd_qt.branch_A1,...
                         bkgd_qt.branch_A4,...
                         bkgd_qt.branch_A5)
    g(:,i)=[ad_d50;ad_d90;ad_udelta;ad_uhat;
          ad_mu;ad_eta;ad_lambda;ad_ahat];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
