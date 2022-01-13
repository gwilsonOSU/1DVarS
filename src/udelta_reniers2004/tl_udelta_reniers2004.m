function [tl_udelta,tl_delta]=tl_udelta_reniers2004(tl_ubar,tl_k,tl_omega,tl_h,tl_Hrms,tl_detady,tl_tau_wind,tl_Dr,tl_fv,tl_ks,tl_d50,bkgd)%,outvar)
%
% TL-code for udelta_reniers2004.m
%
% note, 'bkgd' struct should be taken from 2nd output of
% udelta_reniers2004.m, this contains the NL calculated variables

physicalConstants;

% other fixed params
kappa=0.4;  % von karman
cd=0.002;   % stated in text
fdelta=3;   % stated in text
nintegral=100;  % number of gridpionts for sigma-integration

% % v1, break out NL background vars. This version uses eval() and is a bit slow.
% fld=fields(bkgd);
% for i=1:length(fld)
%   eval([fld{i} ' = bkgd.' fld{i} ';']);
% end

% v2, Break out NL background vars. Note this hard-coded version runs
% significantly faster than if I use eval() to dynamically load all bkgd
% variables
ubar        =bkgd.ubar        ;
k           =bkgd.k           ;
omega       =bkgd.omega       ;
h           =bkgd.h           ;
Hrms        =bkgd.Hrms        ;
detady      =bkgd.detady      ;
Dr          =bkgd.Dr          ;
fv          =bkgd.fv          ;
kabs        =bkgd.kabs        ;
A           =bkgd.A           ;
uorb        =bkgd.uorb        ;
ks          =bkgd.ks          ;
z0          =bkgd.z0          ;
fw          =bkgd.fw          ;
Df          =bkgd.Df          ;
nubar_tb    =bkgd.nubar_tb    ;
Hm0         =bkgd.Hm0         ;
ht          =bkgd.ht          ;
p1          =bkgd.p1          ;
delta       =bkgd.delta       ;
phi_b       =bkgd.phi_b       ;
tau_wind    =bkgd.tau_wind    ;
tau_wind    =bkgd.tau_wind    ;
tau_wave    =bkgd.tau_wave    ;
tau_t       =bkgd.tau_t       ;
nubar_tflow =bkgd.nubar_tflow ;
abs_tau_wind=bkgd.abs_tau_wind;
nubar_twind =bkgd.nubar_twind ;
nubar_twave =bkgd.nubar_twave ;
nubar_t     =bkgd.nubar_t     ;
nu_tsurf    =bkgd.nu_tsurf    ;
sig_s       =bkgd.sig_s       ;
phi_s       =bkgd.phi_s       ;
nums        =bkgd.nums        ;
dens        =bkgd.dens        ;
sig_b       =bkgd.sig_b       ;
sig_0       =bkgd.sig_0       ;
l1d         =bkgd.l1d         ;
l2d         =bkgd.l2d         ;
Ab          =bkgd.Ab          ;
Bbp         =bkgd.Bbp         ;
Cbp         =bkgd.Cbp         ;
Am          =bkgd.Am          ;
alphab_bar  =bkgd.alphab_bar  ;
betab_bar   =bkgd.betab_bar   ;
alpham_bar  =bkgd.alpham_bar  ;
betam_bar   =bkgd.betam_bar   ;
alpha_bar   =bkgd.alpha_bar   ;
beta_bar    =bkgd.beta_bar    ;
F           =bkgd.F           ;
Bb          =bkgd.Bb          ;
Cb          =bkgd.Cb          ;
udelta      =bkgd.udelta      ;

%---------------------------------
% begin TL code
%---------------------------------

% derived params
% kabs=sqrt(k(:,1).^2+k(:,2).^2);
tl_kabs = 1./sqrt(k(:,1).^2+k(:,2).^2).*(k(:,1).*tl_k(:,1)+k(:,2).*tl_k(:,2));
% A=Hrms./(2*sinh(kabs.*h));
tl_A = tl_Hrms./(2*sinh(kabs.*h)) ...
       - Hrms./(2*sinh(kabs.*h)).^2.*( 2*cosh(kabs.*h).*( tl_kabs.*h + kabs.*tl_h ) );
% uorb=A*omega;
tl_uorb = tl_A*omega + A*tl_omega;
% z0=ks/33;
tl_z0 = tl_ks/33;
% fw=1.39*(A./z0).^(-.52);  % eqn (19)
tl_fw = -.52*1.39*(A./z0).^(-.52-1).*( tl_A./z0 - A./z0.^2.*tl_z0 );
% Df=1/(2*sqrt(pi))*rho*fw.*uorb.^3;  % eqn (18)
tl_Df = 1/(2*sqrt(pi))*rho*( tl_fw.*uorb.^3 + 3*fw.*uorb.^2.*tl_uorb );
% nubar_tb=fw.^2.*uorb.^2/(4*omega);  % eqn (21)
tl_nubar_tb = ...
    + 2*fw.*uorb.^2/(4*omega).*tl_fw ...
    + 2*fw.^2.*uorb/(4*omega).*tl_uorb ...
    - fw.^2.*uorb.^2/(4*omega^2)*tl_omega;
% Hm0=Hrms*1.4;
tl_Hm0 = tl_Hrms*1.4;
% ht=h-Hm0/2;  % eqn (3)
tl_ht = tl_h - tl_Hm0/2;
% p1=(A./ks).^(.82);
tl_p1 = .82*(A./ks).^(.82-1)*( tl_A./ks - A./ks.^2.*tl_ks );
% delta=fdelta*0.09*p1.*ks./ht;  % eqn (15)
tl_delta = ...
    + fdelta*0.09*tl_p1.*ks./ht ...
    + fdelta*0.09*p1.*tl_ks./ht ...
    - fdelta*0.09*p1.*ks./ht.^2.*tl_ht;
% phi_b=6./delta.^2;  % eqn (A9)
tl_phi_b = -2*6./delta.^3.*tl_delta;
% tau_wave=Dr.*k/omega;  % eqn (5)
tl_tau_wave = ...
    + tl_Dr.*k/omega ...
    + Dr.*tl_k/omega ...
    - Dr.*k/omega^2*tl_omega;
% tau_t=tau_wind+tau_wave;  % stated in text
tl_tau_t = tl_tau_wind + tl_tau_wave;
% nubar_tflow = 1/6*kappa*ht.*sqrt(g*ht.*abs(detady));  % eqn (11)
if(detady>0)
  sgn=1;
else
  sgn=-1;
end
tl_nubar_tflow = ...
    + 1/6*kappa*tl_ht.*sqrt(g*ht.*abs(detady)) ...
    + .5/6*kappa*ht./sqrt(g*ht.*abs(detady)).*( g*tl_ht.*abs(detady) + sgn*g*ht.*tl_detady );
% abs_tau_wind = sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 );
tl_abs_tau_wind = 1./sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 ).*( tau_wind(:,1).*tl_tau_wind(:,1) + tau_wind(:,2).*tl_tau_wind(:,2) );
% nubar_twind = 1/3*kappa*ht.*sqrt(abs_tau_wind/rho);  % eqn (10)
tl_nubar_twind = ...
    + 1/3*kappa*tl_ht.*sqrt(abs_tau_wind/rho) ...
    + .5/3*kappa*ht./sqrt(abs_tau_wind/rho).*tl_abs_tau_wind/rho;
% nubar_twave = fv*Hrms.*(Dr/rho).^(1/3);  % eqn (9)
tl_nubar_twave = ...
    + tl_fv*Hrms.*(Dr/rho).^(1/3) ...
    + fv*tl_Hrms.*(Dr/rho).^(1/3) ...
    + 1/3*fv*Hrms.*(Dr/rho).^(1/3-1)*tl_Dr/rho;
% nubar_t = sqrt( nubar_tflow.^2 ...
%                 + nubar_twind.^2 ...
%                 + nubar_twave.^2 );  % eqn (12)
tl_nubar_t = 1./sqrt( nubar_tflow.^2 + nubar_twind.^2 + nubar_twave.^2 ).* ...
    ( nubar_tflow.*tl_nubar_tflow + nubar_twind.*tl_nubar_twind + nubar_twave.*tl_nubar_twave );
% nu_tsurf = 3/2*sqrt( nubar_twind.^2 + nubar_twave.^2 );  % eqn (A5)
tl_nu_tsurf = 3/2/sqrt( nubar_twind.^2 + nubar_twave.^2 )*( nubar_twind*tl_nubar_twind + nubar_twave*tl_nubar_twave );
% sig_s = ( nubar_t - 1/3*nu_tsurf ) ...
%         ./( nubar_t - 1/2*nu_tsurf );  % eqn (A4)
tl_sig_s = ...
    + ( tl_nubar_t - 1/3*tl_nu_tsurf )./( nubar_t - 1/2*nu_tsurf ) ...
    - ( nubar_t - 1/3*nu_tsurf )./( nubar_t - 1/2*nu_tsurf )^2*( tl_nubar_t - 1/2*tl_nu_tsurf );
% phi_s = 1./(sig_s/2-1/3);  % eqn (A3)
tl_phi_s = -1./(sig_s/2-1/3).^2.*tl_sig_s/2;
% nums=phi_s.*nubar_t.*sig_s + phi_b.*nubar_tb.*delta;
tl_nums = ...
    + tl_phi_s.*nubar_t.*sig_s ...
    + phi_s.*tl_nubar_t.*sig_s ...
    + phi_s.*nubar_t.*tl_sig_s ...
    + tl_phi_b.*nubar_tb.*delta ...
    + phi_b.*tl_nubar_tb.*delta ...
    + phi_b.*nubar_tb.*tl_delta;
% dens=phi_s.*nubar_t + phi_b.*nubar_tb;
tl_dens = ...
    + tl_phi_s.*nubar_t ...
    + phi_s.*tl_nubar_t ...
    + tl_phi_b.*nubar_tb ...
    + phi_b.*tl_nubar_tb;
% sig_b = nums./dens;  % eqn (A7)
tl_sig_b = tl_nums./dens - nums./dens.^2.*tl_dens;
% sig_0 = z0./ht;  % level at which u=0
tl_sig_0 = tl_z0./ht - z0./ht.^2.*tl_ht;

%--------------------------------------
% Solve for velocity at top of boundary layer, u_delta.  Refer to notes and
% derivations here: ./FgivenUbar_deriv/
%
% Need to do this for each gridpoint separately
%--------------------------------------

% define grid for integration over sigma coordinate
sgridb=linspace(sig_0,delta,nintegral);
sgridbp = linspace(sig_0+tl_sig_0,delta+tl_delta,nintegral);
tl_sgridb = sgridbp - sgridb;
sgridm=linspace(delta,1,nintegral);
sgridmp=linspace(delta+tl_delta,1,nintegral);
tl_sgridm = sgridmp-sgridm;

l1=log(sgridb./sig_0);
tl_l1 = sig_0./sgridb.*( tl_sgridb./sig_0 - sgridb./sig_0.^2.*tl_sig_0 );
l2=log((sig_b-sgridb)./(sig_b-sig_0));
tl_l2 = (sig_b-sig_0)./(sig_b-sgridb).*( (tl_sig_b-tl_sgridb)./(sig_b-sig_0) - (sig_b-sgridb)./(sig_b-sig_0).^2.*(tl_sig_b-tl_sig_0) );
l3=log(sgridm./delta);
tl_l3 = (delta./sgridm).*( tl_sgridm./delta - sgridm./delta.^2.*tl_delta );
l4=log((sig_s-sgridm)./(sig_s-delta));
tl_l4 = (sig_s-delta)./(sig_s-sgridm).*( (tl_sig_s-tl_sgridm)./(sig_s-delta) - (sig_s-sgridm)./(sig_s-delta).^2.*(tl_sig_s-tl_delta) );
l1d=log(delta./sig_0);
tl_l1d = (sig_0./delta).*(tl_delta./sig_0 - delta./sig_0.^2.*tl_sig_0);
l2d=log((sig_b-delta)./(sig_b-sig_0));
tl_l2d = (sig_b-sig_0)./(sig_b-delta).*( (tl_sig_b-tl_delta)./(sig_b-sig_0) - (sig_b-delta)./(sig_b-sig_0).^2.*(tl_sig_b-tl_sig_0) );

for i=1:2  % for each direction

% coef functions for bottom boundary layer velocity profile, see
% ./FgivenUbar_deriv/part1_ubot/
% Ab(i)=ht./(fv*rho.*phi_s.*nubar_t+rho*phi_b.*nubar_tb);  % eqn (B9)
tl_Ab(i) = tl_ht./(fv*rho.*phi_s.*nubar_t+rho*phi_b.*nubar_tb) ...
    - ht./(fv*rho.*phi_s.*nubar_t+rho*phi_b.*nubar_tb).^2.*( ...
        + tl_fv*rho.*phi_s.*nubar_t ...
        + fv*rho.*tl_phi_s.*nubar_t ...
        + fv*rho.*phi_s.*tl_nubar_t ...
        + rho*tl_phi_b.*nubar_tb ...
        + rho*phi_b.*tl_nubar_tb ...
        );
% Bbp(i)=tau_t(i)+Df.*k(i)/omega;
tl_Bbp(i) = tl_tau_t(i) + tl_Df.*k(i)/omega + Df.*tl_k(i)/omega - Df.*k(i)/omega^2*tl_omega;
% Cbp(i)=-Df.*k(i)./(delta*omega);
tl_Cbp(i) = ...
    - tl_Df.*k(i)./(delta*omega) ...
    - Df.*tl_k(i)./(delta*omega) ...
    + Df.*k(i)./(delta*omega)^2*( tl_delta*omega + delta*tl_omega );
alphab=Ab(i).*( Bbp(i)./sig_b.*(l1-l2) - Cbp(i).*l2 );
tl_alphab = ...
    + tl_Ab(i).*( Bbp(i)./sig_b.*(l1-l2) - Cbp(i).*l2 ) ...
    + Ab(i).*( ...
        + tl_Bbp(i)./sig_b.*(l1-l2) ...
        - Bbp(i)./sig_b^2.*(l1-l2)*tl_sig_b ...
        + Bbp(i)./sig_b.*(tl_l1-tl_l2) ...
        - tl_Cbp(i).*l2 ...
        - Cbp(i).*tl_l2 ...
        );
betab=Ab(i).*( -1/sig_b.*(l1-l2) - l2 );
tl_betab = tl_Ab(i).*( -1/sig_b.*(l1-l2) - l2 ) ...
        + Ab(i).*( ...
            + 1/sig_b^2.*(l1-l2)*tl_sig_b ...
            - 1/sig_b.*(tl_l1-tl_l2) ...
            - tl_l2 );

% coef functions for mid-layer velocity profile
% Am(i)=ht./(rho*phi_s.*nubar_t);
tl_Am(i) = ...
    + tl_ht./(rho*phi_s.*nubar_t) ...
    - ht./(rho*phi_s.*nubar_t)^2*(rho*tl_phi_s.*nubar_t + rho*phi_s.*tl_nubar_t);
alpham=alphab(end) + Am(i).*( tau_t(i)./sig_s.*l3 - tau_t(i)./sig_s.*l4 );
tl_alpham = ...
    + tl_alphab(end) ...
    + tl_Am(i).*( tau_t(i)./sig_s.*l3 - tau_t(i)./sig_s.*l4 ) ...
    + Am(i).*( ...
        + tl_tau_t(i)./sig_s.*l3 ...
        - tl_tau_t(i)./sig_s^2.*l3*tl_sig_s ...
        + tau_t(i)./sig_s.*tl_l3 ...
        - tl_tau_t(i)./sig_s.*l4 ...
        + tau_t(i)./sig_s^2.*l4*tl_sig_s ...
        - tau_t(i)./sig_s.*tl_l4 ...
        );
betam=betab(end) + Am(i).*( -1./sig_s.*l3 + (1./sig_s-1).*l4 );
tl_betam = tl_betab(end)...
    + tl_Am(i).*( -1./sig_s.*l3 + (1./sig_s-1).*l4 ) ...
    + Am(i).*( ...
        + 1./sig_s^2.*l3*tl_sig_s ...
        - 1./sig_s.*tl_l3 ...
        - 1./sig_s^2.*l4*tl_sig_s ...
        + (1./sig_s-1).*tl_l4 ...
        );

% define alpha_bar, beta_bar, coefficients for solution ubar.
dsb=diff(sgridb(1:2));
dsm=diff(sgridm(1:2));
tl_dsb = diff(tl_sgridb(1:2));
tl_dsm = diff(tl_sgridm(1:2));
% alphab_bar(i) = trapz(sgridb,alphab);
tl_alphab_bar(i) = trapz(tl_alphab)*dsb + trapz(alphab)*tl_dsb;
% betab_bar(i)  = trapz(sgridb,betab);
tl_betab_bar(i)  = trapz(tl_betab)*dsb + trapz(betab)*tl_dsb;
% alpham_bar(i) = trapz(sgridm,alpham);
tl_alpham_bar(i) = trapz(tl_alpham)*dsm + trapz(alpham)*tl_dsm;
% betam_bar(i)  = trapz(sgridm,betam);
tl_betam_bar(i)  = trapz(tl_betam)*dsm + trapz(betam)*tl_dsm;
% alpha_bar(i)  = alphab_bar(i) + alpham_bar(i);
tl_alpha_bar(i)  = tl_alphab_bar(i) + tl_alpham_bar(i);
% beta_bar(i)   = betab_bar(i) + betam_bar(i);
tl_beta_bar(i)   = tl_betab_bar(i) + tl_betam_bar(i);

% solve for F
% F(i) = (ubar(i) - alpha_bar(i))./beta_bar(i);
tl_F(i) = ...
    + (tl_ubar(i) - tl_alpha_bar(i))./beta_bar(i) ...
    - (ubar(i) - alpha_bar(i))./beta_bar(i).^2.*tl_beta_bar(i);

% evaluate eqns (B9)-(B12) at sigma=delta to get u_delta
% Bb(i) = tau_t(i) - F(i) + Df.*k(i)/omega;  % eqn (B10)
tl_Bb(i) = tl_tau_t(i) - tl_F(i) ...
    + tl_Df.*k(i)/omega ...
    + Df.*tl_k(i)/omega ...
    - Df.*k(i)/omega^2*tl_omega;
% Cb(i)=F(i)-Df.*k(i)./(delta.*omega);  % eqn (B11)
tl_Cb(i) = tl_F(i) ...
    - tl_Df.*k(i)./(delta.*omega) ...
    - Df.*tl_k(i)./(delta.*omega) ...
    + Df.*k(i)./(delta.*omega)^2*(tl_delta.*omega + delta.*tl_omega);

% udelta(i)=Ab(i).*( Bb(i)./sig_b.*l1d - (Bb(i)./sig_b+Cb(i)).*l2d );  % eqn (B12)
tl_udelta(i) = tl_Ab(i).*( Bb(i)./sig_b.*l1d - (Bb(i)./sig_b+Cb(i)).*l2d ) ...
    + Ab(i).*( ...
        + tl_Bb(i)./sig_b.*l1d ...
        - Bb(i)./sig_b^2.*l1d*tl_sig_b ...
        + Bb(i)./sig_b.*tl_l1d ...
        - (tl_Bb(i)./sig_b - Bb(i)./sig_b^2*tl_sig_b + tl_Cb(i)).*l2d ...
        - (Bb(i)./sig_b+Cb(i)).*tl_l2d ...
        );

end  % loop on direction i

% % TEST
% eval(['tl_udelta=tl_' outvar ';']);
