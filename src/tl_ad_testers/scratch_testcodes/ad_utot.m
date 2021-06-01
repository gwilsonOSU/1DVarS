function [ad_h,ad_Hrms,ad_k,ad_udelta]=ad_utot(ad_ut,bkgd)

% break out NL background vars, and recompute some helper variables
omega     =bkgd.omega     ;
t         =bkgd.t         ;
bkgd_uwave=bkgd.bkgd_uwave;
utot      =bkgd.utot      ;
utilde    =bkgd.utilde    ;
udelta    =bkgd.udelta    ;
kabs      =bkgd.kabs      ;
k         =bkgd.k         ;
h         =bkgd.h         ;
Hrms      =bkgd.Hrms      ;
utabs=abs(utilde);
utotabs=abs(utot);
nx=length(h);
nt=length(t);

% physical constants
physicalConstants;
[g,alpha,beta,nu_t]=waveModelParams;

%-----------------------------------
% begin AD code
%-----------------------------------

% init ad vars
ad_udelta=zeros(nx,2);
ad_Uw=zeros(nx,1);
ad_utilde=zeros(nt,nx);
ad_uH=zeros(nt,1);
ad_k=zeros(nx,2);
ad_Hrms=0;

%b03 step through each gridpoint
for i=nx:-1:1

  %8 tl_utot = 1./utot(:,i).*(  ...
  %     utilde(:,i).*tl_utilde ...
  %     + udelta(i,1).*tl_udelta(i,1) ...
  %     + udelta(i,2).*tl_udelta(i,2) );
  % coef=1./utot(:,i);
  % ad_utilde     =ad_utilde     + coef.*utilde(:,i).*ad_ut;
  % ad_udelta(i,1)=ad_udelta(i,1)+ sum(coef.*udelta(i,1).*ad_ut);
  % ad_udelta(i,2)=ad_udelta(i,2)+ sum(coef.*udelta(i,2).*ad_ut);
  % ad_ut=0*ad_ut;
  for j=nt:-1:1
    % tl_utot(j)=1./utot(j,i).*utilde(j,i).*tl_utilde(j,i) ...
    %     + 1./utot(j,i).*udelta(i,1).*tl_udelta(i,1) ...
    %     + 1./utot(j,i).*udelta(i,2).*tl_udelta(i,2);
    ad_utilde(j,i)  =ad_utilde(j,i)  + 1./utot(j,i).*utilde(j,i).*ad_ut(j);
    ad_udelta(i,1)=ad_udelta(i,1)+ 1./utot(j,i).*udelta(i,1).*ad_ut(j);
    ad_udelta(i,2)=ad_udelta(i,2)+ 1./utot(j,i).*udelta(i,2).*ad_ut(j);
    ad_ut(j)=0;
  end

  %2 NOTE: for Hilbert transform, just use brute-force to calculate an
  % equivalent matrix form for imag(hilbert()).  Then manually transpose the
  % matrix to compute the adjoint.
  %
  % See test code in tl_hilbert.m, proving to myself this H matrix works.
  %
  % tl_uH=H*tl_utilde;  % equivalent TL code using the H matrix
  H = imag(hilbert(eye(nt)));  % apply hilbert() to a bunch of delta functions
  ad_utilde(:,i)=ad_utilde(:,i)+H'*ad_uH;
  ad_uH=0*ad_uH;

  %1   tl_utilde=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_h(i),...
  %                                  omega*t,bkgd_uwave(i))';
  [ad_Hmo(i),ad_kabs(i),ad_h(i)]=ad_Uwave_ruessink2012(ad_utilde(:,i)',omega*t,bkgd_uwave(i));
  ad_utilde(:,i)=0*ad_utilde(:,i);

end

%b01 derived params

%3 tl_Uw = omega/2.*( tl_Hrms./sinh(kabs.*h) ...
%                    - Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*(tl_kabs.*h+kabs.*tl_h) );
coef1=omega/2./sinh(kabs.*h);
coef2=omega/2.*Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h);
ad_Hrms=ad_Hrms+ coef1      .*ad_Uw;
ad_kabs=ad_kabs- coef2.*h   .*ad_Uw;
ad_h   =ad_h   - coef2.*kabs.*ad_Uw;
ad_Uw=0;

%2 tl_kabs=1./kabs.*(k(:,1).*tl_k(:,1)+k(:,2).*tl_k(:,2));
coef=1./kabs;
ad_k(:,1)=ad_k(:,1)+coef.*k(:,1).*ad_kabs;  % assume init zero
ad_k(:,2)=ad_k(:,2)+coef.*k(:,2).*ad_kabs;  % assume init zero
ad_kabs=0;
%1 tl_Hmo=tl_Hrms*1.4; % ok
ad_Hrms=ad_Hrms+ad_Hmo*1.4;
