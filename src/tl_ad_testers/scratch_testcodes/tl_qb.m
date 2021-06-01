function tl_q =tl_qb(tl_tanbeta,tl_h,tl_Hrms,tl_kabs,tl_udelta,...
                     tl_Cw,tl_Cc,tl_Cf,bkgd)
%
% TL-code for qtrans_dubarbier.m
%

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

% derived params
tl_Hmo=tl_Hrms*1.4; % ok

% step through each gridpoint
for i=1:nx
  tl_utilde(:,i)=tl_Uwave_ruessink2012(tl_Hmo(i),tl_kabs(i),tl_h(i),...
                                       omega*t,bkgd_uwave(i))';  % ad symmetry ok

  for j=1:length(t)
    tl_utot(j)=1/utot(j,i)*utilde(j,i)*tl_utilde(j,i) ...
        + 1/utot(j,i)*udelta(i,1)*tl_udelta(i,1) ...
        + 1/utot(j,i)*udelta(i,2)*tl_udelta(i,2);
  end
  tl_utot=tl_utot(:);  % ad symmetric
  tl_utabs=sign(utilde(:,i)).*tl_utilde;  % ad symmetric
  tl_utotabs=sign(utot(:,i)).*tl_utot;  % ad symmetric
  tl_qb1=0;
  for j=1:nt
    tl_qb1=tl_qb1+2*utabs(j,i).*tl_utabs(j).*utilde(j,i)/nt ...
                 + utabs(j,i).^2.*tl_utilde(j)/nt;  % ad symmetric
  end
  tl_qb2=0;
  for j=1:nt
    tl_qb2=tl_qb2+2*utabs(j,i).*tl_utabs(j)*udelta(i,1)/nt ...
           + utabs(j,i).^2*tl_udelta(i,1)/nt;  % ad symmetric
  end
  tl_qb3=0;
  for j=1:nt
    tl_qb3=tl_qb3+2*utotabs(j,i).*tl_utotabs(j)/nt;
  end
  tl_q(i) = rho*eps_b/tan_phi*( tl_Cw*qb1(i) ...
                                + Cw*tl_qb1 ...
                                + tl_Cc*qb2(i) ...
                                + Cc*tl_qb2 ...
                                - tl_Cf*tanbeta(i)/tan_phi*qb3(i) ...
                                - Cf*tl_tanbeta(i)/tan_phi*qb3(i) ...
                                - Cf*tanbeta(i)/tan_phi*tl_qb3 );  % NOT symmetric???

end
