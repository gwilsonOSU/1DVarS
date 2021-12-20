%
% TL-AD symmetry test for Uwave_ruessink2012.m
%
clear
addpath ..

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
bkgd.Hmo=waves.H(i)*1.4;
bkgd.k=waves.k(i);
bkgd.omega=waves.sigma;
bkgd.h=waves.h(i);
phs=linspace(0,2*pi,1000);

[Aw,Sw,Uw,bkgd_uwave1]=Uwave_ruessink2012_params(bkgd.Hmo,bkgd.k,bkgd.omega,bkgd.h);
[u,bkgd_uwave2] = Uwave_ruessink2012(phs,Aw,Sw,Uw);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(4,n);  % 1st dim is number of tl input parameters
for i=1:n

  tl_Hmo  =F(1,i);
  tl_k    =F(2,i);
  tl_omega=F(3,i);
  tl_h    =F(4,i);
  % TL model: TL*F
  [tl_Aw,tl_Sw,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_k,tl_omega,tl_h,bkgd_uwave1);
  [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(tl_Aw,tl_Sw,tl_Uw,phs,bkgd_uwave2);
  % AD model: g=AD*(TL*F)
  [ad1_Aw,ad1_Sw,ad1_Uw]=ad_Uwave_ruessink2012(tl_u,tl_r,tl_phi,phs,bkgd_uwave2);
  ad_Aw=tl_Aw+ad1_Aw;
  ad_Sw=tl_Sw+ad1_Sw;
  ad_Uw=tl_Uw+ad1_Uw;
  [ad_Hmo,ad_k,ad_omega,ad_h]=ad_Uwave_ruessink2012_params(ad_Aw,ad_Sw,ad_Uw,bkgd_uwave1);
  g(:,i)=[ad_Hmo; ad_k; ad_omega; ad_h];

end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
