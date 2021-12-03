%
% TL-AD symmetry test for Uwave_ruessink2012.m
%
clear

load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
bkgd.Hmo=waves.H(i)*1.4;
bkgd.k=waves.k(i);
bkgd.omega=waves.sigma;
bkgd.h=waves.h(i);
phs=linspace(0,2*pi,1000);

[u,bkgd_uwave] = Uwave_ruessink2012(phs,bkgd.Hmo,bkgd.k,bkgd.omega,bkgd.h);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(4,n);  % 1st dim is number of tl input parameters
for i=1:n
  % TL model: TL*F
  [tl_u,tl_r,tl_phi]=tl_Uwave_ruessink2012(F(1,i),F(2,i),F(3,i),F(4,i),phs,bkgd_uwave);
  % AD model: g=AD*(TL*F)
  [ad_Hmo,ad_k,ad_omega,ad_h]=ad_Uwave_ruessink2012(tl_u,tl_r,tl_phi,phs,bkgd_uwave);
  g(:,i)=[ad_Hmo; ad_k; ad_omega; ad_h];
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
