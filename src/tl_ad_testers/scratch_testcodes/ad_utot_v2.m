function [ad_h,ad_Hrms,ad_k,ad_udelta]=ad_utot_v2(ad_ut,bkgd)
%
% TEST-CODE
%

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

ad_utilde=zeros(nt,nx);
ad_Uw=0;
ad_udelta=zeros(nx,2);
ad_uH=zeros(nt,1);
ad_Hrms=0;
ad_k=zeros(nx,2);

% step through each gridpoint
for i=nx:-1:1

  for j=nt:-1:1
    %3   tl_ut(j)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
    %       + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
    %       + 1/utot(j,i)*udelta(i,2)*tl_udelta(i,2);
    ad_utilde(j,i)=ad_utilde(j,i)+1/utot(j,i)*utilde(j,i)*ad_ut(j);
    ad_udelta(i,1)=ad_udelta(i,1)+1/utot(j,i)*udelta(i,1)*ad_ut(j);
    ad_udelta(i,2)=ad_udelta(i,2)+1/utot(j,i)*udelta(i,2)*ad_ut(j);
    ad_ut(j)=0;
  end

  H = imag(hilbert(eye(nt)));  % apply hilbert() to a bunch of delta functions
  %2 tl_uH=H*tl_utilde(:,i);
  ad_utilde(:,i)=ad_utilde(:,i)+H'*ad_uH;
  ad_uH=0*ad_uH;

  %1 tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_h(i),...
  %                                      omega*t,bkgd_uwave(i))';  % ad symmetry ok
  [ad_Hmo(i),ad_kabs(i),ad_h(i)]=ad_Uwave_ruessink2012(ad_utilde(:,i)',omega*t,bkgd_uwave(i));
  
end


%3 tl_Uw = omega/2.*( tl_Hrms./sinh(kabs.*h) ...
%                    - Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*(tl_kabs.*h+kabs.*tl_h) );
ad_Hrms=ad_Hrms+ omega/2./sinh(kabs.*h).*ad_Uw;
ad_kabs=ad_kabs- omega/2.*Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*h   .*ad_Uw;
ad_h   =ad_h   - omega/2.*Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*kabs.*ad_Uw;
ad_Uw=0;

%2 tl_kabs=1./kabs.*(k(:,1).*tl_k(:,1)+k(:,2).*tl_k(:,2));
ad_k(:,1)=ad_k(:,1)+ 1./kabs.*k(:,1).*ad_kabs;
ad_k(:,2)=ad_k(:,2)+ 1./kabs.*k(:,2).*ad_kabs;
ad_kabs=0;

%1 tl_Hmo=tl_Hrms*1.4; % ok
ad_Hrms=ad_Hrms+ad_Hmo*1.4;
ad_Hmo=0;
