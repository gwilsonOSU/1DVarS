%
% TL-AD symmetry test for Uwave_ruessink2012.m
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
% inoutvar='utabs';
% inoutvar='utotabs';
% inoutvar='qb1';
% % inoutvar='qb1';
% % inoutvar='qb2';
% % inoutvar='qb2';
% % inoutvar='qb3';
% % inoutvar='qb3';
% inoutvar='qb';
% % inoutvar='qs1';
% % inoutvar='qs2';
% % inoutvar='qs3';
% inoutvar='qs';
% inoutvar='q';
% % inoutvar='Q';

load ../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
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
bkgd.d50=200e-6;
bkgd.kabs=waves.k(i);
bkgd.k=waves.k(i)*[cos(waves.theta(i)) sin(waves.theta(i))];
bkgd.udelta=udelta_reniers2004(bkgd.ubar,bkgd.k,bkgd.omega,bkgd.h,...
                               bkgd.Hrms,bkgd.detady,bkgd.windW,...
                               bkgd.Dr,bkgd.fv,bkgd.d50);
bkgd.ws=ws_brownLawler(bkgd.d50);
bkgd.Cw=0.00483;
bkgd.Cc=0.02002;
bkgd.Cf=0.01173;
bkgd.Ka=0.631e-4;
bkgd.tanbeta=0.01;
if(length(i)==1)
  bkgd.x=bkgd.tanbeta;  % wonky code hack
end

% bkgd state
[q,qb,qs,qa,bkgd_qt]=qtrans_dubarbier(bkgd.tanbeta,...
                                      bkgd.h      ,...
                                      bkgd.Hrms   ,...
                                      bkgd.kabs   ,...
                                      bkgd.omega  ,...
                                      bkgd.udelta ,...
                                      bkgd.ws     ,...
                                      bkgd.Cw     ,...
                                      bkgd.Cc     ,...
                                      bkgd.Cf     ,...
                                      bkgd.Ka     );

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F = eps*rand(11,n);  % 1st dim is number of tl input parameters
for i=1:n
  % TL model: TL*F
  tl_tanbeta=F(1,i);
  tl_h      =F(2,i);
  tl_Hrms   =F(3,i);
  tl_kabs   =F(4,i);
  tl_udelta =F(5:6,i)';
  tl_ws     =F(7,i);
  tl_Cw     =F(8,i);
  tl_Cc     =F(9,i);
  tl_Cf     =F(10,i);
  tl_Ka     =F(11,i);
  tl_Q=tl_qtrans_dubarbier(tl_tanbeta,...
                           tl_h,...
                           tl_Hrms,...
                           tl_kabs,...
                           tl_udelta,...
                           tl_ws,...
                           tl_Cw,...
                           tl_Cc,...
                           tl_Cf,...
                           tl_Ka,...
                           bkgd_qt);
  % AD model: g=AD*(TL*F)
  [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_udelta,ad_ws,...
   ad_Cw,ad_Cc,ad_Cf,ad_Ka] = ...
      ad_qtrans_dubarbier(tl_Q,bkgd_qt);
  g(:,i)=[ad_tanbeta(:);ad_h(:);ad_Hrms(:);ad_kabs(:);...
          ad_udelta(:);ad_ws(:);...
          ad_Cw(:);ad_Cc(:);ad_Cf(:);ad_Ka];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
