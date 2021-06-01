function tl_ut=tl_utot(tl_h,tl_Hrms,tl_k,tl_udelta,bkgd)
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

% derived params
tl_Hmo=tl_Hrms*1.4; % ok
tl_kabs=1./kabs.*(k(:,1).*tl_k(:,1)+k(:,2).*tl_k(:,2));
tl_Uw = omega/2.*( tl_Hrms./sinh(kabs.*h) ...
                   - Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*(tl_kabs.*h+kabs.*tl_h) );

% step through each gridpoint
for i=1:nx

  % intra-wave velocity information, using Ruessink et al. (2012).  Note,
  % comparing to Hsu et al. (2006) it seems as though Dubarbier et al. are
  % using \tilde{U} and \tilde{u} (upper and lower case) interchangably for
  % "wave velocity".  There are a few other obvious sloppy notational errors
  % like using "omega_s" for fall velocity
  tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_h(i),...
                                       omega*t,bkgd_uwave(i))';  % ad symmetry ok

  H = imag(hilbert(eye(nt)));  % apply hilbert() to a bunch of delta functions
  % tl_uH=imag(hilbert(tl_utilde));
  tl_uH=H*tl_utilde(:,i);  % ad symmetry ok
  for j=1:nt
    tl_ut(j)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
        + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
        + 1/utot(j,i)*udelta(i,2)*tl_udelta(i,2);
  end

end

