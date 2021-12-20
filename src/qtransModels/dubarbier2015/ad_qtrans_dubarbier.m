function [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,...
                     ad_Cw,ad_Cc,ad_Cf,ad_Ka] = ...
    ad_qtrans_dubarbier(ad_Q,bkgd)%,invar)
%
% AD-code for tl_qtrans_dubarbier.m
%

% break out NL background vars, and recompute some helper variables
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
utabs=abs(utilde);
utotabs=abs(utot);
nx=length(h);
nt=length(t);

% physical constants
physicalConstants;

%-----------------------------------
% begin AD code
%-----------------------------------

% init ad vars
ad_q=zeros(nx,1);
ad_qs=zeros(nx,1);
ad_qb=zeros(nx,1);
ad_qa=zeros(nx,1);
ad_Cw=0;
ad_qs1=0;
ad_qs2=0;
ad_qs3=0;
ad_qb1=0;
ad_qb2=0;
ad_qb3=0;
ad_Cc=0;
ad_Cf=0;
ad_tanbeta=zeros(nx,1);
ad_ws=zeros(nx,1);
ad_udelta=zeros(nx,2);
ad_Ka=0;
ad_Au=zeros(nx,1);
ad_Aw=0;
ad_Uw=zeros(nx,1);
ad_utilde=zeros(nt,nx);
ad_utot=zeros(nt,nx);
ad_utotabs=zeros(nt,nx);
ad_utabs=zeros(nt,nx);
ad_uH=zeros(nt,1);
ad_Hmo=zeros(nx,1);
ad_Hrms=zeros(nx,1);
ad_h=zeros(nx,1);
ad_kabs=zeros(nx,1);
ad_Au_nums=0;
ad_Au_dens=zeros(nx,1);
ad_Au_dens_1=zeros(nx,1);
ad_omega=0;

% % TEST
% eval(['ad_' invar '=ad_Q;']);
% if(~strcmp(invar,'Q'))
%   ad_Q=zeros(nx,1);
% end

%b04 normalize to get volumetric transport, m2/s units (bed volume per unit
% time per unit width), and convert to row vectors

%4 tl_Q =tl_q(:)/qnorm;
ad_q =ad_q+ad_Q/qnorm;
ad_Q=0;

%b03 step through each gridpoint
for i=nx:-1:1

  %b03-2 energetics model, eqns (11)-(14)

  %12 tl_q(i)=tl_qa(i)+tl_qb(i)+tl_qs(i);
  ad_qa(i)=ad_qa(i)+ad_q(i);
  ad_qb(i)=ad_qb(i)+ad_q(i);
  ad_qs(i)=ad_qs(i)+ad_q(i);
  ad_q(i)=0;

  %11 tl_qs(i) = rho*eps_s/ws(i)*( tl_Cw*qs1(i) ...
  %                        + Cw*tl_qs1 ...
  %                        + tl_Cc*qs2(i) ...
  %                        + Cc*tl_qs2 ...
  %                        - tl_Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) ...
  %                        - Cf*eps_s*tl_tanbeta(i)/ws(i)*qs3(i) ...
  %                        + Cf*eps_s*tanbeta(i)/ws(i)^2*qs3(i)*tl_ws(i) ...
  %                        - Cf*eps_s*tanbeta(i)/ws(i)*tl_qs3 ) ...
  %         - rho*eps_s/ws(i)^2*( Cw*qs1(i) ...
  %                            + Cc*qs2(i) ...
  %                            - Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) )*tl_ws(i);
  coef1=rho*eps_s/ws(i);
  coef2=coef1/ws(i)*( Cw*qs1(i) + Cc*qs2(i) - Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) );
  ad_Cw        =ad_Cw        + coef1*qs1(i)                         *ad_qs(i);
  ad_qs1       =ad_qs1       + coef1*Cw                             *ad_qs(i);
  ad_Cc        =ad_Cc        + coef1*qs2(i)                         *ad_qs(i);
  ad_qs2       =ad_qs2       + coef1*Cc                             *ad_qs(i);
  ad_Cf        =ad_Cf        - coef1*eps_s*tanbeta(i)/ws(i)*qs3(i)     *ad_qs(i);
  ad_tanbeta(i)=ad_tanbeta(i)- coef1*Cf*eps_s/ws(i)*qs3(i)             *ad_qs(i);
  ad_ws(i)        =ad_ws(i)        + coef1*Cf*eps_s*tanbeta(i)/ws(i)^2*qs3(i)*ad_qs(i);
  ad_qs3       =ad_qs3       - coef1*Cf*eps_s*tanbeta(i)/ws(i)         *ad_qs(i);
  ad_ws(i)        =ad_ws(i)        - coef2                               .*ad_qs(i);
  ad_qs(i)=0;

  %10 tl_qs3=mean( 5*utotabs(:,i).^4.*tl_utotabs );
  ad_utotabs(:,i)=ad_utotabs(:,i)+5*utotabs(:,i).^4.*ad_qs3/nt;
  ad_qs3=0;

  %9 tl_qs2=mean( 3*utotabs(:,i).^2.*tl_utotabs*udelta(i,1) ...
  %              + utotabs(:,i).^3*tl_udelta(i,1) );
  ad_utotabs(:,i)=ad_utotabs(:,i)+3*utotabs(:,i).^2*udelta(i,1)/nt*ad_qs2;
  ad_udelta(i,1)=ad_udelta(i,1)+mean(utotabs(:,i).^3)*ad_qs2;
  ad_qs2=0;

  %8 tl_qs1=mean( 3*utabs(:,i).^2.*tl_utabs.*utilde(:,i) ...
  %              + utabs(:,i).^3.*tl_utilde);
  ad_utabs(:,i)=ad_utabs(:,i)+3*utabs(:,i).^2.*utilde(:,i)*ad_qs1/nt;
  ad_utilde(:,i)=ad_utilde(:,i)+utabs(:,i).^3*ad_qs1/nt;
  ad_qs1=0;

  %9 tl_qb(i) = rho*eps_b/tan_phi*( tl_Cw*qb1(i) ...
  %                               + Cw*tl_qb1 ...
  %                               + tl_Cc*qb2(i) ...
  %                               + Cc*tl_qb2 ...
  %                               - tl_Cf*tanbeta(i)/tan_phi*qb3(i) ...
  %                               - Cf*tl_tanbeta(i)/tan_phi*qb3(i) ...
  %                               - Cf*tanbeta(i)/tan_phi*tl_qb3 );  % NOT symmetric???
  coef=rho*eps_b/tan_phi;
  ad_Cw =ad_Cw + coef*qb1(i)*ad_qb(i);
  ad_qb1=ad_qb1+ coef*Cw    *ad_qb(i);
  ad_Cc =ad_Cc + coef*qb2(i)*ad_qb(i);
  ad_qb2=ad_qb2+ coef*Cc    *ad_qb(i);
  ad_Cf        =ad_Cf        - coef*tanbeta(i)/tan_phi*qb3(i)*ad_qb(i);
  ad_tanbeta(i)=ad_tanbeta(i)- coef*Cf/tan_phi*qb3(i)        *ad_qb(i);
  ad_qb3       =ad_qb3       - coef*Cf*tanbeta(i)/tan_phi    *ad_qb(i);
  ad_qb(i)=0;

  for j=nt:-1:1
    %6   tl_qb3=tl_qb3+3*utotabs(j,i)^2.*tl_utotabs(j)/nt;
    ad_utotabs(j,i)=ad_utotabs(j,i)+ 3*utotabs(j,i)^2/nt*ad_qb3;
  end
  ad_qb3=0;

  for j=nt:-1:1
    %5   tl_qb2=tl_qb2+2*utabs(j,i).*tl_utabs(j)*udelta(i,1)/nt ...
    %          + utabs(j,i).^2*tl_udelta(i,1)/nt;  % ad symmetric
    ad_utabs(j,i)   =ad_utabs(j,i)   + 2*utabs(j,i)*udelta(i,1)/nt*ad_qb2;
    ad_udelta(i,1)=ad_udelta(i,1)+ utabs(j,i).^2/nt           *ad_qb2;
  end
  ad_qb2=0;

  for j=nt:-1:1
    %4   tl_qb1=tl_qb1+2*utabs(j,i).*tl_utabs(j).*utilde(j,i)/nt ...
    %                + utabs(j,i).^2.*tl_utilde(j)/nt;  % ad symmetric
    ad_utabs(j,i) =ad_utabs(j,i) + 2*utabs(j,i).*utilde(j,i)/nt*ad_qb1;
    ad_utilde(j,i)=ad_utilde(j,i)+ utabs(j,i).^2/nt            *ad_qb1;
  end
  ad_qb1=0;

  %3 tl_utotabs=sign(utot(:,i)).*tl_utot;  % ad symmetric
  ad_utot(:,i)=ad_utot(:,i)+sign(utot(:,i)).*ad_utotabs(:,i);
  ad_utotabs(:,i)=0;

  %2 tl_utabs(:,i)=sign(utilde(:,i)).*tl_utilde(:,i);  % ad symmetric
  ad_utilde(:,i)=ad_utilde(:,i)+sign(utilde(:,i)).*ad_utabs(:,i);
  ad_utabs(:,i)=0;

  %1 tl_qa(i) = -tl_Ka*Au(i)*Aw(i) ...
  %         -Ka*tl_Au*Aw(i) ...
  %         -Ka*Au(i)*tl_Aw;
  ad_Ka=ad_Ka- Au(i)*Aw(i)*ad_qa(i);
  ad_Au(i)=ad_Au(i)- Ka*Aw(i)   *ad_qa(i);
  ad_Aw=ad_Aw- Ka*Au(i)   *ad_qa(i);
  ad_qa(i)=0;

  %b03-1 intra-wave velocity information, using Ruessink et al. (2012).  Note,
  % comparing to Hsu et al. (2006) it seems as though Dubarbier et al. are
  % using \tilde{U} and \tilde{u} (upper and lower case) interchangably for
  % "wave velocity".  There are a few other obvious sloppy notational errors
  % like using "omega_s" for fall velocity

  for j=nt:-1:1
    % tl_utot(j,i)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
    %     + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
    %     + 1/utot(j,i)*udelta(i,2)*tl_udelta(i,2);
    ad_utilde(j,i)=ad_utilde(j,i)+ 1/utot(j,i)*utilde(j,i)*ad_utot(j,i);
    ad_udelta(i,1)=ad_udelta(i,1)+ 1/utot(j,i)*udelta(i,1)*ad_utot(j,i);
    ad_udelta(i,2)=ad_udelta(i,2)+ 1/utot(j,i)*udelta(i,2)*ad_utot(j,i);
    ad_utot(j,i)=0;
  end

  %7 tl_Aw = omega*tl_Uw(i) + tl_omega*Uw(i);
  ad_Uw(i)=ad_Uw(i)+ omega*ad_Aw;
  ad_omega = ad_omega + Uw(i)*ad_Aw;
  ad_Aw=0;

  %6 tl_Au = tl_Au_nums./Au_dens(i) - Au_nums(i)./Au_dens(i).^2.*tl_Au_dens;
  ad_Au_nums   =ad_Au_nums   + 1/Au_dens(i)           *ad_Au(i);
  ad_Au_dens(i)=ad_Au_dens(i)- Au_nums(i)/Au_dens(i)^2*ad_Au(i);
  ad_Au(i)=0;

  %5 tl_Au_dens = (3/2)*Au_dens_1(i).^(3/2-1).*tl_Au_dens_1;
  ad_Au_dens_1(i)=ad_Au_dens_1(i)+(3/2)*Au_dens_1(i).^(3/2-1).*ad_Au_dens(i);
  ad_Au_dens(i)=0;

  %4 tl_Au_dens_1 = mean( 2*uH(:,i).*tl_uH );
  ad_uH=ad_uH+2*uH(:,i).*ad_Au_dens_1(i)/nt;
  ad_Au_dens_1(i)=0;

  %3 tl_Au_nums = mean( 3*uH(:,i).^2.*tl_uH );
  % ad_uH=ad_uH+3*uH(:,i).^2.*ad_Au_nums/nt;
  for j=nt:-1:1
    % tl_Au_nums = tl_Au_nums + 3*uH(:,i).^2.*tl_uH/nt;
    ad_uH(j)=ad_uH(j)+3*uH(j,i)^2/nt*ad_Au_nums;
  end
  ad_Au_nums=0;

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

  %1 tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_omega,tl_h(i),...
  %                                      omega*t,bkgd_uwave(i))';
  [ad1_Hmo,ad1_kabs,ad1_omega,ad1_h]=ad_Uwave_ruessink2012(ad_utilde(:,i)',0,0,omega*t,0,0,bkgd_uwave(i));
  ad_Hmo(i) =ad_Hmo(i) +ad1_Hmo;
  ad_kabs(i)=ad_kabs(i)+ad1_kabs;
  ad_omega = ad_omega + ad1_omega;
  ad_h(i)   =ad_h(i)   +ad1_h;
  ad_utilde(:,i)=0;

end

%b01 derived params

%3 tl_Uw = omega/2.*( tl_Hrms./sinh(kabs.*h) ...
%                    - Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h).*(tl_kabs.*h+kabs.*tl_h) ) ...
%         + tl_omega.*Hrms./(2*sinh(kabs.*h));
coef1=omega/2./sinh(kabs.*h);
coef2=omega/2.*Hrms./sinh(kabs.*h).^2.*cosh(kabs.*h);
ad_Hrms=ad_Hrms+ coef1      .*ad_Uw;
ad_kabs=ad_kabs- coef2.*h   .*ad_Uw;
ad_h   =ad_h   - coef2.*kabs.*ad_Uw;
ad_omega = ad_omega + sum(Hrms./(2*sinh(kabs.*h)).*ad_Uw);
ad_Uw=0;

% tl_Hmo=tl_Hrms*1.4; % ok
ad_Hrms=ad_Hrms+1.4*ad_Hmo;
ad_Hmo=0;
