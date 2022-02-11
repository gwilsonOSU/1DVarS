function tl_Q =tl_qtrans_dubarbier(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_ws,tl_Aw,tl_Sw,tl_Uw,...
                     tl_Cw,tl_Cc,tl_Cf,tl_Ka,bkgd)%,outvar)
%
% TL-code for qtrans_dubarbier.m
%

% break out NL background vars, and recompute some helper variables.  Note,
% is much faster to hard-code this rather than use 'fields(bkgd)' to
% auto-populate a list.
h          =bkgd.h          ;
Hrms       =bkgd.Hrms       ;
omega      =bkgd.omega      ;
udelta     =bkgd.udelta     ;
ws         =bkgd.ws         ;
Cw         =bkgd.Cw         ;
Cc         =bkgd.Cc         ;
Cf         =bkgd.Cf         ;
Ka         =bkgd.Ka         ;
eps_b      =bkgd.eps_b      ;
eps_s      =bkgd.eps_s      ;
tan_phi    =bkgd.tan_phi    ;
nt         =bkgd.nt         ;
nx         =bkgd.nx         ;
Hmo        =bkgd.Hmo        ;
kabs       =bkgd.kabs       ;
Uwq        =bkgd.Uwq        ;
tanbeta    =bkgd.tanbeta    ;
t          =bkgd.t          ;
utilde     =bkgd.utilde     ;
bkgd_uwave =bkgd.bkgd_uwave ;
uH         =bkgd.uH         ;
Au_nums    =bkgd.Au_nums    ;
Au_dens_1  =bkgd.Au_dens_1  ;
Au_dens    =bkgd.Au_dens    ;
Au         =bkgd.Au         ;
Awq        =bkgd.Awq        ;
utot       =bkgd.utot       ;
qb1        =bkgd.qb1        ;
qb2        =bkgd.qb2        ;
qb3        =bkgd.qb3        ;
qs1        =bkgd.qs1        ;
qs2        =bkgd.qs2        ;
qs3        =bkgd.qs3        ;
qa         =bkgd.qa         ;
qs         =bkgd.qs         ;
qb         =bkgd.qb         ;
qnorm      =bkgd.qnorm      ;
Qa         =bkgd.Qa         ;
Qs         =bkgd.Qs         ;
Qb         =bkgd.Qb         ;
Aw         =bkgd.Aw         ;
Sw         =bkgd.Sw         ;
Uw         =bkgd.Uw         ;

% define some stray NL vars not saved in bkgd struct
utabs=abs(utilde);
utotabs=abs(utot);
nt=length(t);
nx=length(h);

% physical constants
physicalConstants;

% derived params
tl_Hmo=tl_Hrms*1.4; % ok
tl_Uwq= omega/2.*( tl_Hrms./sinh(kabs.*h) ...
                   - Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*(tl_kabs.*h+kabs.*tl_h) ) ...
        + tl_omega.*Hrms./(2*sinh(kabs.*h));

% step through each gridpoint
for i=1:nx

  % for masked points, make a dummy output with all fields set to 0
  if(Hrms(i)==0)
  tl_q(i)  =0;
  else

  % intra-wave velocity information, using Ruessink et al. (2012).  Note,
  % comparing to Hsu et al. (2006) it seems as though Dubarbier et al. are
  % using \tilde{U} and \tilde{u} (upper and lower case) interchangably for
  % "wave velocity".  There are a few other obvious sloppy notational errors
  % like using "omega_s" for fall velocity
  tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Aw(i),tl_Sw(i),tl_Uw(i),...
                                       omega*t,bkgd_uwave(i))';

  H = imag(hilbert(eye(length(t))));  % apply hilbert() to a bunch of delta functions
  % tl_uH=imag(hilbert(tl_utilde));
  tl_uH=H*tl_utilde(:,i);
  tl_Au_nums = mean( 3*uH(:,i).^2.*tl_uH );
  tl_Au_dens_1(i) = mean( 2*uH(:,i).*tl_uH );
  tl_Au_dens(i) = (3/2)*Au_dens_1(i).^(3/2-1).*tl_Au_dens_1(i);
  tl_Au(i) = tl_Au_nums./Au_dens(i) - Au_nums(i)./Au_dens(i).^2.*tl_Au_dens(i);
  tl_Awq(i) = omega*tl_Uwq(i) + tl_omega*Uwq(i);
  for j=1:length(t)
    tl_utot(j,i)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
        + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
        + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,2);
  end

  % energetics model, eqns (11)-(14)
  tl_qa(i) = -tl_Ka*Au(i)*Awq(i) ...
          -Ka*tl_Au(i)*Awq(i) ...
          -Ka*Au(i)*tl_Awq(i);  % ad symmetric
  tl_utabs(:,i)=sign(utilde(:,i)).*tl_utilde(:,i);  % ad symmetric
  tl_utotabs(:,i)=sign(utot(:,i)).*tl_utot(:,i);  % ad symmetric

  tl_qb1=0;
  for j=1:nt
    tl_qb1=tl_qb1+2*utabs(j,i).*tl_utabs(j,i).*utilde(j,i)/nt ...
                 + utabs(j,i).^2.*tl_utilde(j,i)/nt;  % ad symmetric
  end
  tl_qb2=0;
  for j=1:nt
    tl_qb2=tl_qb2+2*utabs(j,i).*tl_utabs(j,i)*udelta(i,1)/nt ...
           + utabs(j,i).^2*tl_udelta(i,1)/nt;  % ad symmetric
  end
  tl_qb3=0;
  for j=1:nt
    tl_qb3=tl_qb3+3*utotabs(j,i).^2.*tl_utotabs(j,i)/nt;
  end
  tl_qb(i) = rho*eps_b/tan_phi*( tl_Cw*qb1(i) ...
                              + Cw*tl_qb1 ...
                              + tl_Cc*qb2(i) ...
                              + Cc*tl_qb2 ...
                              - tl_Cf*tanbeta(i)/tan_phi*qb3(i) ...
                              - Cf*tl_tanbeta(i)/tan_phi*qb3(i) ...
                              - Cf*tanbeta(i)/tan_phi*tl_qb3 );  % NOT symmetric???
  tl_qs1=mean( 3*utabs(:,i).^2.*tl_utabs(:,i).*utilde(:,i) ...
               + utabs(:,i).^3.*tl_utilde(:,i));  % ad symmetric
  tl_qs2=mean( 3*utotabs(:,i).^2.*tl_utotabs(:,i)*udelta(i,1) ...
               + utotabs(:,i).^3*tl_udelta(i,1) );  % ad symmetric
  tl_qs3=mean( 5*utotabs(:,i).^4.*tl_utotabs(:,i) );  % ad symmetric
  tl_qs(i) = rho*eps_s/ws(i)*( tl_Cw*qs1(i) ...
                         + Cw*tl_qs1 ...
                         + tl_Cc*qs2(i) ...
                         + Cc*tl_qs2 ...
                         - tl_Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) ...
                         - Cf*eps_s*tl_tanbeta(i)/ws(i)*qs3(i) ...
                         + Cf*eps_s*tanbeta(i)/ws(i)^2*qs3(i)*tl_ws(i) ...
                         - Cf*eps_s*tanbeta(i)/ws(i)*tl_qs3 ) ...
          - rho*eps_s/ws(i)^2*( Cw*qs1(i) ...
                             + Cc*qs2(i) ...
                             - Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) )*tl_ws(i);  % ad symmetric
  tl_q(i)=tl_qa(i)+tl_qb(i)+tl_qs(i);

  end  % masking

end

% normalize to get volumetric transport, m2/s units (bed volume per unit
% time per unit width), and convert to row vectors
tl_Q =tl_q(:)/qnorm;

% % TEST
% eval(['tl_Q=tl_' outvar ';']);
