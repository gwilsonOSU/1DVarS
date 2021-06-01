%
% TL-AD symmetry test for ws_brownLawler.m
%
clear

d50=200e-6;

ws=ws_brownLawler(d50);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=10;
F = eps*rand(1,n);  % 1st dim is number of tl input parameters
for i=1:n
  % TL model: TL*F
  tl_ws=tl_ws_brownLawler(F(i),d50);
  % AD model: g=AD*(TL*F)
  g(i)=ad_ws_brownLawler(tl_ws,d50);
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
