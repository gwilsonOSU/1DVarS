%
% TL-AD symmetry test for qb.m
%
addpath src_jtech2018/waveModel/
clear

load src_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
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

% bkgd state
[q,bkgd_qt]=qb(bkgd.tanbeta,...
               bkgd.h      ,...
               bkgd.Hrms   ,...
               bkgd.kabs   ,...
               bkgd.omega  ,...
               bkgd.udelta ,...
               bkgd.Cw     ,...
               bkgd.Cc     ,...
               bkgd.Cf     );

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F = eps*rand(9,n);  % 1st dim is number of tl input parameters
% F(2,:)=0;
% F(3,:)=0;
% F(4,:)=0;
for i=1:n
  % TL model: TL*F
  tl_tanbeta=F(1,i);
  tl_h      =F(2,i);
  tl_Hrms   =F(3,i);
  tl_kabs   =F(4,i);
  tl_udelta =F(5:6,i)';
  tl_Cw     =F(7,i);
  tl_Cc     =F(8,i);
  tl_Cf     =F(9,i);
  tl_q=tl_qb(tl_tanbeta,...
             tl_h,...
             tl_Hrms,...
             tl_kabs,...
             tl_udelta,...
             tl_Cw,...
             tl_Cc,...
             tl_Cf,...
             bkgd_qt);
  % AD model: g=AD*(TL*F)
  [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_udelta,ad_Cw,ad_Cc,ad_Cf] =ad_qb(tl_q,bkgd_qt);
  g(:,i)=[ad_tanbeta(:);ad_h(:);ad_Hrms(:);ad_kabs(:);
          ad_udelta(:);ad_Cw(:);ad_Cc(:);ad_Cf(:)];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
