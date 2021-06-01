function [ad_tanbeta,ad_h,ad_Hrms,ad_kabs,ad_udelta,ad_Cw,ad_Cc,ad_Cf] =ad_qb(ad_q,bkgd)


% break out NL background vars, and recompute some helper variables
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end
utabs=abs(utilde);
utotabs=abs(utot);
nt=length(t);
nx=length(h);

% physical constants
physicalConstants;

%-----------------------------------

ad_Cw=0;
ad_Cc=0;
ad_Cf=0;
ad_udelta=[0 0];
ad_kabs=0;
ad_Hrms=0;
ad_h=0;
ad_tanbeta=0;
ad_qb1=0;
ad_qb2=0;
ad_qb3=0;
ad_utotabs=zeros(nt,1);
ad_utot=zeros(nt,1);
ad_utabs=zeros(nt,1);
ad_utilde=zeros(nt,1);

% step through each gridpoint
for i=1:nx

  %9 tl_q(i) = rho*eps_b/tan_phi*( tl_Cw*qb1(i) ...
  %                               + Cw*tl_qb1 ...
  %                               + tl_Cc*qb2(i) ...
  %                               + Cc*tl_qb2 ...
  %                               - tl_Cf*tanbeta(i)/tan_phi*qb3(i) ...
  %                               - Cf*tl_tanbeta(i)/tan_phi*qb3(i) ...
  %                               - Cf*tanbeta(i)/tan_phi*tl_qb3 );  % NOT symmetric???
  coef=rho*eps_b/tan_phi;
  ad_Cw =ad_Cw + coef*qb1(i)*ad_q(i);
  ad_qb1=ad_qb1+ coef*Cw    *ad_q(i);
  ad_Cc =ad_Cc + coef*qb2(i)*ad_q(i);
  ad_qb2=ad_qb2+ coef*Cc    *ad_q(i);
  ad_Cf        =ad_Cf        - coef*tanbeta(i)/tan_phi*qb3(i)*ad_q(i);
  ad_tanbeta(i)=ad_tanbeta(i)- coef*Cf/tan_phi*qb3(i)        *ad_q(i);
  ad_qb3       =ad_qb3       - coef*Cf*tanbeta(i)/tan_phi    *ad_q(i);
  ad_q(i)=0;

  for j=nt:-1:1
    %8   tl_qb3=tl_qb3+2*utotabs(j,i).*tl_utotabs(j)/nt;
    ad_utotabs(j)=ad_utotabs(j)+ 2*utotabs(j,i)/nt*ad_qb3;
  end
  ad_qb3=0;

  for j=nt:-1:1
    %7   tl_qb2=tl_qb2+2*utabs(j,i).*tl_utabs(j)*udelta(i,1)/nt ...
    %          + utabs(j,i).^2*tl_udelta(i,1)/nt;  % ad symmetric
    ad_utabs(j)   =ad_utabs(j)   + 2*utabs(j,i)*udelta(i,1)/nt*ad_qb2;
    ad_udelta(i,1)=ad_udelta(i,1)+ utabs(j,i).^2/nt           *ad_qb2;
  end
  ad_qb2=0;

  for j=nt:-1:1
    %6   tl_qb1=tl_qb1+2*utabs(j,i).*tl_utabs(j).*utilde(j,i)/nt ...
    %                + utabs(j,i).^2.*tl_utilde(j)/nt;  % ad symmetric
    ad_utabs(j) =ad_utabs(j) + 2*utabs(j,i).*utilde(j,i)/nt*ad_qb1;
    ad_utilde(j)=ad_utilde(j)+ utabs(j,i).^2/nt            *ad_qb1;
  end
  ad_qb1=0;

  %5 tl_utotabs=sign(utot(:,i)).*tl_utot;  % ad symmetric
  ad_utot=ad_utot+sign(utot(:,i)).*ad_utotabs;
  ad_utotabs=0;

  %4 tl_utabs=sign(utilde(:,i)).*tl_utilde;  % ad symmetric
  ad_utilde=ad_utilde+sign(utilde(:,i)).*ad_utabs;
  ad_utabs=0;

  for j=nt:-1:1
    %2   tl_utot(j)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
    %       + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
    %       + 1/utot(j,i)*udelta(i,2)*tl_udelta(i,2);
    ad_utilde(j,i)=ad_utilde(j,i)+ 1/utot(j,i)*utilde(j,i)*ad_utot(j);
    ad_udelta(i,1)=ad_udelta(i,1)+ 1/utot(j,i)*udelta(i,1)*ad_utot(j);
    ad_udelta(i,2)=ad_udelta(i,2)+ 1/utot(j,i)*udelta(i,2)*ad_utot(j);
    ad_utot(j)=0;
  end

  %1 tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_h(i),...
  %                                      omega*t,bkgd_uwave(i))';  % ad symmetry ok
  % 
  [ad_Hmo(i),ad_kabs(i),ad_h(i)]=ad_Uwave_ruessink2012(ad_utilde(:,i)',omega*t,bkgd_uwave(i));
  ad_utilde(:,i)=0;

end

% tl_Hmo=tl_Hrms*1.4; % ok
ad_Hrms=ad_Hrms+1.4*ad_Hmo;
ad_Hmo=0;

