function [ad_ubar,ad_k,ad_omega,ad_h,ad_Hrms,ad_detady,ad_tau_wind,ad_Dr,ad_fv,ad_ks,ad_d50]=ad_udelta_reniers2004(ad_udelta,ad_delta,bkgd)%,invar)
%
% AD-code for tl_udelta_reniers2004.m
%

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

% init ad vars
ad_Ab=zeros(1,2);
ad_Bb=zeros(1,2);
ad_Cb=zeros(1,2);
ad_Bbp=zeros(1,2);
ad_Cbp=zeros(1,2);
ad_Am=zeros(1,2);
ad_sig_b=0;
ad_l1d=0;
ad_l2d=0;
ad_sig_0=0;
ad_F=zeros(1,2);
ad_Df=0;
ad_k=zeros(1,2);
ad_tau_t=zeros(1,2);
ad_ubar=zeros(1,2);
ad_alpha_bar=zeros(1,2);
ad_beta_bar=zeros(1,2);
ad_alphab_bar=zeros(1,2);
ad_betab_bar=zeros(1,2);
ad_alpham_bar=zeros(1,2);
ad_betam_bar=zeros(1,2);
ad_sig_s=0;
ad_ht=0;
ad_phi_s=0;
ad_nubar_t=0;
ad_fv=0;
ad_phi_b=0;
ad_nubar_tb=0;
ad_z0=0;
ad_nums=0;
ad_dens=0;
ad_nu_tsurf=0;
ad_nubar_twind=0;
ad_nubar_twave=0;
ad_nubar_tflow=0;
ad_Hrms=0;
ad_Dr=0;
ad_abs_tau_wind=0;
ad_tau_wind=zeros(1,2);
ad_windW=zeros(1,2);
ad_detady=0;
ad_tau_wave=0;
ad_w1=0;
ad_windW=zeros(1,2);
ad_phi=0;
ad_p1=0;
ad_ks=0;
ad_A=0;
ad_h=0;
ad_Hm0=0;
ad_fw=0;
ad_uorb=0;
ad_d50=0;
ad_kabs=0;
ad_l1=zeros(1,nintegral);
ad_l2=zeros(1,nintegral);
ad_l3=zeros(1,nintegral);
ad_l4=zeros(1,nintegral);
ad_omega=0;
ad_dsm=0;
ad_dsb=0;
ad_betam=zeros(1,nintegral);
ad_alpham=zeros(1,nintegral);
ad_betab=zeros(1,nintegral);
ad_alphab=zeros(1,nintegral);
ad_sgridm=zeros(1,nintegral);
ad_sgridb=zeros(1,nintegral);

% % TEST
% eval(['ad_' invar '=ad_udelta;']);
% if(~strcmp(invar,'udelta'))
%   eval(['ad_udelta=zeros(1,2);']);
% end

for i=2:-1:1  % for each direction

% NL-vars needed for this block
dsb = (delta-sig_0)/(nintegral-1);
dsm = (1-delta)/(nintegral-1);
for n=1:nintegral
  sgridb(n) = sig_0 + (n-1)*dsb;
  sgridm(n) = delta + (n-1)*dsm;
end
l1=log(sgridb./sig_0);
l2=log((sig_b-sgridb)./(sig_b-sig_0));
l3=log(sgridm./delta);
l4=log((sig_s-sgridm)./(sig_s-delta));
l1d=log(delta./sig_0);
l2d=log((sig_b-delta)./(sig_b-sig_0));
alphab=Ab(i)*( Bbp(i)/sig_b*(l1-l2) - Cbp(i)*l2 );
betab=Ab(i)*( -1/sig_b*(l1-l2) - l2 );
alpham=alphab(end) + Am(i).*( tau_t(i)/sig_s.*l3 - tau_t(i)/sig_s.*l4 );
betam=betab(end) + Am(i).*( -1/sig_s.*l3 + (1/sig_s-1).*l4 );
if(detady>0)
  sgn=1;
else
  sgn=-1;
end

% evaluate eqns (B9)-(B12) at sigma=delta to get u_delta
%53 tl_udelta(i) = tl_Ab(i).*( Bb(i)./sig_b.*l1d - (Bb(i)./sig_b+Cb(i)).*l2d ) ...
%     + Ab(i).*( ...
%         + tl_Bb(i)./sig_b.*l1d ...
%         - Bb(i)./sig_b^2.*l1d*tl_sig_b ...
%         + Bb(i)./sig_b.*tl_l1d ...
%         - (tl_Bb(i)./sig_b - Bb(i)./sig_b^2*tl_sig_b + tl_Cb(i)).*l2d ...
%         - (Bb(i)./sig_b+Cb(i)).*tl_l2d ...
%         );
ad_Ab(i)=ad_Ab(i)+ ( Bb(i)./sig_b.*l1d - (Bb(i)./sig_b+Cb(i)).*l2d ).*ad_udelta(i);
ad_Bb(i)=ad_Bb(i)+ Ab(i)./sig_b.*l1d                                .*ad_udelta(i);
ad_sig_b=ad_sig_b- Ab(i).*Bb(i)./sig_b^2.*l1d                       .*ad_udelta(i);
ad_l1d  =ad_l1d  + Ab(i).*Bb(i)./sig_b                              .*ad_udelta(i);
ad_Bb(i)=ad_Bb(i)- Ab(i).*l2d./sig_b                                .*ad_udelta(i);
ad_sig_b=ad_sig_b+ Ab(i).*l2d.*Bb(i)./sig_b^2                       .*ad_udelta(i);
ad_Cb(i)=ad_Cb(i)- Ab(i).*l2d                                       .*ad_udelta(i);
ad_l2d  =ad_l2d  - Ab(i).*(Bb(i)./sig_b+Cb(i))                      .*ad_udelta(i);
ad_udelta(i)=0;
%52 tl_Cb(i) = tl_F(i) ...
%     - tl_Df.*k(i)./(delta.*omega) ...
%     - Df.*tl_k(i)./(delta.*omega) ...
%     + Df.*k(i)./(delta.*omega)^2*(tl_delta.*omega + delta.*tl_omega);
ad_F(i) =ad_F(i) + 1                              *ad_Cb(i);
ad_Df   =ad_Df   - k(i)/(delta.*omega)            *ad_Cb(i);
ad_k(i) =ad_k(i) - Df/(delta.*omega)              *ad_Cb(i);
ad_delta=ad_delta+ Df*k(i)./(delta.*omega)^2*omega*ad_Cb(i);
ad_omega=ad_omega+ Df*k(i)./(delta.*omega)^2*delta*ad_Cb(i);
ad_Cb(i)=0;
%51 tl_Bb(i) = tl_tau_t(i) - tl_F(i) ...
%     + tl_Df.*k(i)/omega ...
%     + Df.*tl_k(i)/omega ...
%     - Df.*k(i)/omega^2*tl_omega;
ad_tau_t(i)=ad_tau_t(i)+ 1              *ad_Bb(i);
ad_F(i)    =ad_F(i)    - 1              *ad_Bb(i);
ad_Df      =ad_Df      + k(i)/omega     *ad_Bb(i);
ad_k(i)    =ad_k(i)    + Df/omega       *ad_Bb(i);
ad_omega   =ad_omega   - Df*k(i)/omega^2*ad_Bb(i);
ad_Bb(i)=0;

% solve for F
%50 tl_F(i) = ...
%     + (tl_ubar(i) - tl_alpha_bar(i))./beta_bar(i) ...
%     - (ubar(i) - alpha_bar(i))./beta_bar(i).^2.*tl_beta_bar(i);
ad_ubar(i)     =ad_ubar(i)     + 1/beta_bar(i)                         *ad_F(i);
ad_alpha_bar(i)=ad_alpha_bar(i)- 1/beta_bar(i)                         *ad_F(i);
ad_beta_bar(i) =ad_beta_bar(i) - (ubar(i) - alpha_bar(i))/beta_bar(i)^2*ad_F(i);
ad_F(i)=0;

% define alpha_bar, beta_bar, coefficients for solution ubar.
%49 tl_beta_bar(i)   = tl_betab_bar(i) + tl_betam_bar(i);
ad_betab_bar(i)=ad_betab_bar(i)+ ad_beta_bar(i);
ad_betam_bar(i)=ad_betam_bar(i)+ ad_beta_bar(i);
ad_beta_bar(i)=0;
%48 tl_alpha_bar(i)  = tl_alphab_bar(i) + tl_alpham_bar(i);
ad_alphab_bar(i)=ad_alphab_bar(i)+ ad_alpha_bar(i);
ad_alpham_bar(i)=ad_alpham_bar(i)+ ad_alpha_bar(i);
ad_alpha_bar(i)=0;

%47d tl_betam_bar(i) = tl_betam_bar(i) + tl_dsm/2*sum(betam(2:end));
%47c tl_betam_bar(i) = tl_betam_bar(i) + tl_dsm/2*sum(betam(1:end-1));
%47b tl_betam_bar(i) = tl_betam_bar(i) + dsm/2*sum(tl_betam(2:end));
%47a tl_betam_bar(i) = dsm/2*sum(tl_betam(1:end-1));
ad_dsm=ad_dsm+ 1/2*sum(betam(2:end))*ad_betam_bar(i);
ad_dsm=ad_dsm+1/2*sum(betam(1:end-1))*ad_betam_bar(i);
ad_betam(2:nintegral) = ad_betam(2:nintegral) + dsm/2*ad_betam_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_betam(1:nintegral-1) = ad_betam(1:nintegral-1) + dsm/2*ad_betam_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_betam_bar(i)=0;

%46d tl_alpham_bar(i) = tl_alpham_bar(i) + tl_dsm/2*sum(alpham(2:end));
%46c tl_alpham_bar(i) = tl_alpham_bar(i) + tl_dsm/2*sum(alpham(1:end-1));
%46b tl_alpham_bar(i) = tl_alpham_bar(i) + dsm/2*sum(tl_alpham(2:end));
%46a tl_alpham_bar(i) = dsm/2*sum(tl_alpham(1:end-1));
ad_dsm=ad_dsm+ 1/2*sum(alpham(2:end))*ad_alpham_bar(i);
ad_dsm=ad_dsm+1/2*sum(alpham(1:end-1))*ad_alpham_bar(i);
ad_alpham(2:nintegral) = ad_alpham(2:nintegral) + dsm/2*ad_alpham_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_alpham(1:nintegral-1) = ad_alpham(1:nintegral-1) + dsm/2*ad_alpham_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_alpham_bar(i)=0;

%45d tl_betab_bar(i) = tl_betab_bar(i) + tl_dsb/2*sum(betab(2:end));
%45c tl_betab_bar(i) = tl_betab_bar(i) + tl_dsb/2*sum(betab(1:end-1));
%45b tl_betab_bar(i) = tl_betab_bar(i) + dsb/2*sum(tl_betab(2:end));
%45a tl_betab_bar(i) = dsb/2*sum(tl_betab(1:end-1));
ad_dsb=ad_dsb+ 1/2*sum(betab(2:end))*ad_betab_bar(i);
ad_dsb=ad_dsb+1/2*sum(betab(1:end-1))*ad_betab_bar(i);
ad_betab(2:nintegral) = ad_betab(2:nintegral) + dsb/2*ad_betab_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_betab(1:nintegral-1) = ad_betab(1:nintegral-1) + dsb/2*ad_betab_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_betab_bar(i)=0;

%44d tl_alphab_bar(i) = tl_alphab_bar(i) + tl_dsb/2*sum(alphab(2:end));
%44c tl_alphab_bar(i) = tl_alphab_bar(i) + tl_dsb/2*sum(alphab(1:end-1));
%44b tl_alphab_bar(i) = tl_alphab_bar(i) + dsb/2*sum(tl_alphab(2:end));
%44a tl_alphab_bar(i) = dsb/2*sum(tl_alphab(1:end-1));
ad_dsb=ad_dsb+ 1/2*sum(alphab(2:end))*ad_alphab_bar(i);
ad_dsb=ad_dsb+1/2*sum(alphab(1:end-1))*ad_alphab_bar(i);
ad_alphab(2:nintegral) = ad_alphab(2:nintegral) + dsb/2*ad_alphab_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_alphab(1:nintegral-1) = ad_alphab(1:nintegral-1) + dsb/2*ad_alphab_bar(i)*ones(1,nintegral-1);  % distribute sum over array
ad_alphab_bar(i)=0;

% coef functions for mid-layer velocity profile
%43 tl_betam = tl_betab(end)...
%     + tl_Am(i).*( -1/sig_s.*l3 + (1/sig_s-1).*l4 ) ...
%     + Am(i).*( ...
%         + 1/sig_s^2.*l3*tl_sig_s ...
%         - 1/sig_s.*tl_l3 ...
%         - 1/sig_s^2.*l4*tl_sig_s ...
%         + (1/sig_s-1).*tl_l4 ...
%         );
ad_betab(end)=ad_betab(end)+ sum(1                               .*ad_betam);
ad_Am(i)     =ad_Am(i)     + sum(( -1/sig_s*l3 + (1/sig_s-1)*l4 ).*ad_betam);
ad_sig_s     =ad_sig_s     + sum(Am(i)*1/sig_s^2*l3              .*ad_betam);
ad_l3        =ad_l3        - Am(i)*1/sig_s                   .*ad_betam;
ad_sig_s     =ad_sig_s     - sum(Am(i)*1/sig_s^2*l4              .*ad_betam);
ad_l4        =ad_l4        + Am(i)*(1/sig_s-1)               .*ad_betam;
ad_betam=zeros(1,nintegral);
%42 tl_alpham = ...
%     + tl_alphab(end) ...
%     + tl_Am(i)*( tau_t(i)/sig_s*l3 - tau_t(i)/sig_s*l4 ) ...
%     + Am(i)*( ...
%         + tl_tau_t(i)/sig_s*l3 ...
%         - tau_t(i)/sig_s^2*l3*tl_sig_s ...
%         + tau_t(i)/sig_s*tl_l3 ...
%         - tl_tau_t(i)/sig_s*l4 ...
%         + tau_t(i)/sig_s^2*l4*tl_sig_s ...
%         - tau_t(i)/sig_s*tl_l4 ...
%         );
ad_alphab(end)=ad_alphab(end)+ sum(1                                        *ad_alpham);
ad_Am(i)      =ad_Am(i)      + sum(( tau_t(i)/sig_s*l3 - tau_t(i)/sig_s*l4 ).*ad_alpham);
ad_tau_t(i)   =ad_tau_t(i)   + sum(Am(i)/sig_s*l3                           .*ad_alpham);
ad_sig_s      =ad_sig_s      - sum(Am(i)*tau_t(i)/sig_s^2*l3                .*ad_alpham);
ad_l3         =ad_l3         + Am(i)*tau_t(i)/sig_s                     .*ad_alpham;
ad_tau_t(i)   =ad_tau_t(i)   - sum(Am(i)/sig_s*l4                           .*ad_alpham);
ad_sig_s      =ad_sig_s      + sum(Am(i)*tau_t(i)/sig_s^2*l4                .*ad_alpham);
ad_l4         =ad_l4         - Am(i)*tau_t(i)/sig_s                     .*ad_alpham;
ad_alpham=zeros(1,nintegral);
%41 tl_Am(i) = ...
%     + tl_ht/(rho*phi_s*nubar_t) ...
%     - ht/(rho*phi_s*nubar_t)^2*(rho*tl_phi_s*nubar_t + rho*phi_s*tl_nubar_t);
ad_ht     =ad_ht     + 1/(rho*phi_s*nubar_t)               *ad_Am(i);
ad_phi_s  =ad_phi_s  - ht/(rho*phi_s*nubar_t)^2*rho*nubar_t*ad_Am(i);
ad_nubar_t=ad_nubar_t- ht/(rho*phi_s*nubar_t)^2*rho*phi_s  *ad_Am(i);
ad_Am(i)=0;
%40 tl_betab = tl_Ab(i)*( -1/sig_b*(l1-l2) - l2 ) ...
%         + Ab(i)*( ...
%             + 1/sig_b^2*(l1-l2)*tl_sig_b ...
%             - 1/sig_b*(tl_l1-tl_l2) ...
%             - tl_l2 );
ad_Ab(i)=ad_Ab(i)+ sum(( -1/sig_b*(l1-l2) - l2 ).*ad_betab);
ad_sig_b=ad_sig_b+ sum(Ab(i)*1/sig_b^2*(l1-l2)  .*ad_betab);
ad_l1   =ad_l1   - Ab(i)*1/sig_b            .*ad_betab;
ad_l2   =ad_l2   + Ab(i)*1/sig_b            .*ad_betab;
ad_l2   =ad_l2   - Ab(i)                    .*ad_betab;
ad_betab=zeros(1,nintegral);
%39 tl_alphab = ...
%     + tl_Ab(i)*( Bbp(i)/sig_b*(l1-l2) - Cbp(i)*l2 ) ...
%     + Ab(i)*( ...
%         + tl_Bbp(i)/sig_b*(l1-l2) ...
%         - Bbp(i)/sig_b^2*(l1-l2)*tl_sig_b ...
%         + Bbp(i)/sig_b*(tl_l1-tl_l2) ...
%         - tl_Cbp(i)*l2 ...
%         - Cbp(i)*tl_l2 ...
%         );
ad_Ab(i) =ad_Ab(i) + sum(( Bbp(i)/sig_b*(l1-l2) - Cbp(i)*l2 ).*ad_alphab);
ad_Bbp(i)=ad_Bbp(i)    + sum(Ab(i)/sig_b*(l1-l2)             .*ad_alphab);
ad_sig_b =ad_sig_b     - sum(Ab(i)*Bbp(i)/sig_b^2*(l1-l2)    .*ad_alphab);
ad_l1    =ad_l1        + (Ab(i)*Bbp(i)/sig_b              .*ad_alphab);
ad_l2    =ad_l2        - (Ab(i)*Bbp(i)/sig_b              .*ad_alphab);
ad_Cbp(i)=ad_Cbp(i)    - sum(Ab(i)*l2                        .*ad_alphab);
ad_l2    =ad_l2        - (Ab(i)*Cbp(i)                    .*ad_alphab);
ad_alphab=zeros(1,nintegral);
%38 tl_Cbp(i) = ...
%     - tl_Df*k(i)/(delta*omega) ...
%     - Df*tl_k(i)/(delta*omega) ...
%     + Df*k(i)/(delta*omega)^2*( tl_delta*omega + delta*tl_omega );
ad_Df   =ad_Df   - k(i)/(delta*omega)           *ad_Cbp(i);
ad_k(i) =ad_k(i) - Df/(delta*omega)             *ad_Cbp(i);
ad_delta=ad_delta+ Df*k(i)/(delta*omega)^2*omega*ad_Cbp(i);
ad_omega=ad_omega+ Df*k(i)/(delta*omega)^2*delta*ad_Cbp(i);
ad_Cbp(i)=0;
%37 tl_Bbp(i) = tl_tau_t(i) + tl_Df*k(i)/omega + Df*tl_k(i)/omega - Df*k(i)/omega^2*tl_omega;
ad_tau_t(i)=ad_tau_t(i)+ 1              *ad_Bbp(i);
ad_Df      =ad_Df      + k(i)/omega     *ad_Bbp(i);
ad_k(i)    =ad_k(i)    + Df/omega       *ad_Bbp(i);
ad_omega   =ad_omega   - Df*k(i)/omega^2*ad_Bbp(i);
ad_Bbp(i)=0;
%36 tl_Ab(i) = tl_ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb) ...
%     - ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*( ...
%         + tl_fv*rho*phi_s*nubar_t ...
%         + fv*rho*tl_phi_s*nubar_t ...
%         + fv*rho*phi_s*tl_nubar_t ...
%         + rho*tl_phi_b*nubar_tb ...
%         + rho*phi_b*tl_nubar_tb ...
%         );
ad_ht      =ad_ht      + 1/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb)                      *ad_Ab(i);
ad_fv      =ad_fv      - ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*rho*phi_s*nubar_t*ad_Ab(i);
ad_phi_s   =ad_phi_s   - ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*fv*rho*nubar_t   *ad_Ab(i);
ad_nubar_t =ad_nubar_t - ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*fv*rho*phi_s     *ad_Ab(i);
ad_phi_b   =ad_phi_b   - ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*rho*nubar_tb     *ad_Ab(i);
ad_nubar_tb=ad_nubar_tb- ht/(fv*rho*phi_s*nubar_t+rho*phi_b*nubar_tb).^2*rho*phi_b        *ad_Ab(i);
ad_Ab(i)=0;

end  % loop on direction i

% grid
%35 tl_l2d = (sig_b-sig_0)./(sig_b-delta).*( (tl_sig_b-tl_delta)./(sig_b-sig_0) - (sig_b-delta)./(sig_b-sig_0).^2.*(tl_sig_b-tl_sig_0) );
ad_sig_b=ad_sig_b+ (sig_b-sig_0)./(sig_b-delta)./(sig_b-sig_0)                  .*ad_l2d;
ad_delta=ad_delta- (sig_b-sig_0)./(sig_b-delta)./(sig_b-sig_0)                  .*ad_l2d;
ad_sig_b=ad_sig_b- (sig_b-sig_0)./(sig_b-delta).*(sig_b-delta)./(sig_b-sig_0).^2.*ad_l2d;
ad_sig_0=ad_sig_0+ (sig_b-sig_0)./(sig_b-delta).*(sig_b-delta)./(sig_b-sig_0).^2.*ad_l2d;
ad_l2d=0;
%34 tl_l1d = (sig_0./delta).*(tl_delta./sig_0 - delta./sig_0.^2.*tl_sig_0);
ad_delta=ad_delta+ (sig_0./delta)./sig_0          .*ad_l1d;
ad_sig_0=ad_sig_0- (sig_0./delta).*delta./sig_0.^2.*ad_l1d;
ad_l1d=0;
%33 tl_l4 = (sig_s-delta)./(sig_s-sgridm).*( (tl_sig_s-tl_sgridm)./(sig_s-delta) - (sig_s-sgridm)./(sig_s-delta).^2.*(tl_sig_s-tl_delta) );
ad_sig_s =ad_sig_s + sum((sig_s-delta)./(sig_s-sgridm)./(sig_s-delta)                   .*ad_l4);
ad_sgridm=ad_sgridm- ((sig_s-delta)./(sig_s-sgridm)./(sig_s-delta)                   .*ad_l4);
ad_sig_s =ad_sig_s - sum((sig_s-delta)./(sig_s-sgridm).*(sig_s-sgridm)./(sig_s-delta).^2.*ad_l4);
ad_delta =ad_delta + sum((sig_s-delta)./(sig_s-sgridm).*(sig_s-sgridm)./(sig_s-delta).^2.*ad_l4);
ad_l4=zeros(1,nintegral);
%32 tl_l3 = (delta./sgridm).*( tl_sgridm./delta - sgridm./delta.^2.*tl_delta );
ad_sgridm=ad_sgridm+ (delta./sgridm)./delta           .*ad_l3;
ad_delta =ad_delta - sum((delta./sgridm).*sgridm./delta.^2.*ad_l3);
ad_l3=zeros(1,nintegral);
%31 tl_l2 = (sig_b-sig_0)./(sig_b-sgridb).*( (tl_sig_b-tl_sgridb)./(sig_b-sig_0) - (sig_b-sgridb)./(sig_b-sig_0).^2.*(tl_sig_b-tl_sig_0) );
ad_sig_b =ad_sig_b + sum((sig_b-sig_0)./(sig_b-sgridb)./(sig_b-sig_0)                   .*ad_l2);
ad_sgridb=ad_sgridb- ((sig_b-sig_0)./(sig_b-sgridb)./(sig_b-sig_0)                   .*ad_l2);
ad_sig_b =ad_sig_b - sum((sig_b-sig_0)./(sig_b-sgridb).*(sig_b-sgridb)./(sig_b-sig_0).^2.*ad_l2);
ad_sig_0 =ad_sig_0 + sum((sig_b-sig_0)./(sig_b-sgridb).*(sig_b-sgridb)./(sig_b-sig_0).^2.*ad_l2);
ad_l2=zeros(1,nintegral);
%30 tl_l1 = sig_0./sgridb.*( tl_sgridb./sig_0 - sgridb./sig_0.^2.*tl_sig_0 );
ad_sgridb=ad_sgridb+ sig_0./sgridb./sig_0           .*ad_l1;
ad_sig_0 =ad_sig_0 - sum(sig_0./sgridb.*sgridb./sig_0.^2.*ad_l1);
ad_l1=zeros(1,nintegral);
for n=1:nintegral
  %29a2   tl_sgridm(n) = tl_delta + (n-1)*tl_dsm;
  ad_delta=ad_delta+ 1    *ad_sgridm(n);
  ad_dsm  =ad_dsm  + (n-1)*ad_sgridm(n);
  ad_sgridm(n)=0;
  %29a1   tl_sgridb(n) = tl_sig_0 + (n-1)*tl_dsb;
  ad_sig_0=ad_sig_0+ 1    *ad_sgridb(n);
  ad_dsb  =ad_dsb  + (n-1)*ad_sgridb(n);
  ad_sgridb(n)=0;
end
%28 tl_dsm = -tl_delta/(nintegral-1);
ad_delta=ad_delta- 1/(nintegral-1)*ad_dsm;
ad_dsm=0;
%27 tl_dsb = (tl_delta-tl_sig_0)/(nintegral-1);
ad_delta=ad_delta+ 1/(nintegral-1)*ad_dsb;
ad_sig_0=ad_sig_0- 1/(nintegral-1)*ad_dsb;
ad_dsb=0;

% derived params
%26 tl_sig_0 = tl_z0./ht - z0./ht.^2.*tl_ht;
ad_z0=ad_z0+ 1./ht    .*ad_sig_0;
ad_ht=ad_ht- z0./ht.^2.*ad_sig_0;
ad_sig_0=0;
%25 tl_sig_b = tl_nums./dens - nums./dens.^2.*tl_dens;
ad_nums=ad_nums+ 1./dens      .*ad_sig_b;
ad_dens=ad_dens- nums./dens.^2.*ad_sig_b;
ad_sig_b=0;
%24 tl_dens = ...
%     + tl_phi_s.*nubar_t ...
%     + phi_s.*tl_nubar_t ...
%     + tl_phi_b.*nubar_tb ...
%     + phi_b.*tl_nubar_tb;
ad_phi_s   =ad_phi_s   + nubar_t .*ad_dens;
ad_nubar_t =ad_nubar_t + phi_s   .*ad_dens;
ad_phi_b   =ad_phi_b   + nubar_tb.*ad_dens;
ad_nubar_tb=ad_nubar_tb+ phi_b   .*ad_dens;
ad_dens=0;
%23 tl_nums = ...
%     + tl_phi_s.*nubar_t.*sig_s ...
%     + phi_s.*tl_nubar_t.*sig_s ...
%     + phi_s.*nubar_t.*tl_sig_s ...
%     + tl_phi_b.*nubar_tb.*delta ...
%     + phi_b.*tl_nubar_tb.*delta ...
%     + phi_b.*nubar_tb.*tl_delta;
ad_phi_s   =ad_phi_s   + nubar_t.*sig_s .*ad_nums;
ad_nubar_t =ad_nubar_t + phi_s.*sig_s   .*ad_nums;
ad_sig_s   =ad_sig_s   + phi_s.*nubar_t .*ad_nums;
ad_phi_b   =ad_phi_b   + nubar_tb.*delta.*ad_nums;
ad_nubar_tb=ad_nubar_tb+ phi_b.*delta   .*ad_nums;
ad_delta   =ad_delta   + phi_b.*nubar_tb.*ad_nums;
ad_nums=0;
%22 tl_phi_s = -1./(sig_s/2-1/3).^2.*tl_sig_s/2;
ad_sig_s=ad_sig_s- 1./(sig_s/2-1/3).^2/2.*ad_phi_s;
ad_phi_s=0;
%21 tl_sig_s = ...
%     + ( tl_nubar_t - 1/3*tl_nu_tsurf )./( nubar_t - 1/2*nu_tsurf ) ...
%     - ( nubar_t - 1/3*nu_tsurf )./( nubar_t - 1/2*nu_tsurf )^2*( tl_nubar_t - 1/2*tl_nu_tsurf );
ad_nubar_t =ad_nubar_t + 1./( nubar_t - 1/2*nu_tsurf )    .*ad_sig_s;
ad_nu_tsurf=ad_nu_tsurf- 1./( nubar_t - 1/2*nu_tsurf ).*1/3*ad_sig_s;
ad_nubar_t =ad_nubar_t - ( nubar_t - 1/3*nu_tsurf )./( nubar_t - 1/2*nu_tsurf )^2    *ad_sig_s;
ad_nu_tsurf=ad_nu_tsurf+ ( nubar_t - 1/3*nu_tsurf )./( nubar_t - 1/2*nu_tsurf )^2*1/2*ad_sig_s;
ad_sig_s=0;
%20 tl_nu_tsurf = 3/2/sqrt( nubar_twind.^2 + nubar_twave.^2 )*( nubar_twind*tl_nubar_twind + nubar_twave*tl_nubar_twave );
ad_nubar_twind=ad_nubar_twind+ 3/2/sqrt( nubar_twind.^2 + nubar_twave.^2 )*nubar_twind*ad_nu_tsurf;
ad_nubar_twave=ad_nubar_twave+ 3/2/sqrt( nubar_twind.^2 + nubar_twave.^2 )*nubar_twave*ad_nu_tsurf;
ad_nu_tsurf=0;
%19 tl_nubar_t = 1./sqrt( nubar_tflow.^2 + nubar_twind.^2 + nubar_twave.^2 ).* ...
%     ( nubar_tflow.*tl_nubar_tflow + nubar_twind.*tl_nubar_twind + nubar_twave.*tl_nubar_twave );
ad_nubar_tflow=ad_nubar_tflow+ 1./sqrt( nubar_tflow.^2 + nubar_twind.^2 + nubar_twave.^2 ).*nubar_tflow.*ad_nubar_t;
ad_nubar_twind=ad_nubar_twind+ 1./sqrt( nubar_tflow.^2 + nubar_twind.^2 + nubar_twave.^2 ).*nubar_twind.*ad_nubar_t;
ad_nubar_twave=ad_nubar_twave+ 1./sqrt( nubar_tflow.^2 + nubar_twind.^2 + nubar_twave.^2 ).*nubar_twave.*ad_nubar_t;
ad_nubar_t=0;
%18 tl_nubar_twave = ...
%     + tl_fv*Hrms.*(Dr/rho).^(1/3) ...
%     + fv*tl_Hrms.*(Dr/rho).^(1/3) ...
%     + 1/3*fv*Hrms.*(Dr/rho).^(1/3-1)*tl_Dr/rho;
ad_fv  =ad_fv  + Hrms.*(Dr/rho).^(1/3)             *ad_nubar_twave;
ad_Hrms=ad_Hrms+ fv.*(Dr/rho).^(1/3)               *ad_nubar_twave;
ad_Dr  =ad_Dr  + 1/3*fv*Hrms.*(Dr/rho).^(1/3-1)/rho*ad_nubar_twave;
ad_nubar_twave=0;
%17 tl_nubar_twind = ...
%     + 1/3*kappa*tl_ht.*sqrt(abs_tau_wind/rho) ...
%     + .5/3*kappa*ht./sqrt(abs_tau_wind/rho).*tl_abs_tau_wind/rho;
ad_ht          =ad_ht          + 1/3*kappa.*sqrt(abs_tau_wind/rho)        *ad_nubar_twind;
ad_abs_tau_wind=ad_abs_tau_wind+ .5/3*kappa*ht./sqrt(abs_tau_wind/rho)/rho*ad_nubar_twind;
ad_nubar_twind=0;
%16 tl_abs_tau_wind = 1./sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 ).*( tau_wind(:,1).*tl_tau_wind(:,1) + tau_wind(:,2).*tl_tau_wind(:,2) );
ad_tau_wind(:,1)=ad_tau_wind(:,1)+ 1./sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 ).*tau_wind(:,1).*ad_abs_tau_wind;
ad_tau_wind(:,2)=ad_tau_wind(:,2)+ 1./sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 ).*tau_wind(:,2).*ad_abs_tau_wind;
ad_abs_tau_wind=0;
%15 tl_nubar_tflow = ...
%     + 1/6*kappa*tl_ht.*sqrt(g*ht.*abs(detady)) ...
%     + .5/6*kappa*ht./sqrt(g*ht.*abs(detady)).*( g*tl_ht.*abs(detady) + sgn*g*ht.*tl_detady );
ad_ht    =ad_ht    + 1/6*kappa.*sqrt(g*ht.*abs(detady))                    *ad_nubar_tflow;
ad_ht    =ad_ht    + .5/6*kappa*ht./sqrt(g*ht.*abs(detady)).*g.*abs(detady)*ad_nubar_tflow;
ad_detady=ad_detady+ .5/6*kappa*ht./sqrt(g*ht.*abs(detady)).*sgn*g*ht     .*ad_nubar_tflow;
ad_nubar_tflow=0;
%14 tl_tau_t = tl_tau_wind + tl_tau_wave;
ad_tau_wind=ad_tau_wind+ ad_tau_t;
ad_tau_wave=ad_tau_wave+ ad_tau_t;
ad_tau_t=0;
%13 tl_tau_wave = ...
%     + tl_Dr.*k/omega ...
%     + Dr.*tl_k/omega ...
%     - Dr.*k/omega^2*tl_omega;
ad_Dr   =ad_Dr   + sum(k/omega      .*ad_tau_wave);
ad_k    =ad_k    + Dr/omega     *ad_tau_wave;
ad_omega=ad_omega- sum(Dr.*k/omega^2.*ad_tau_wave);
ad_tau_wave=0;
%12 tl_phi_b = -2*6./delta.^3.*tl_delta;
ad_delta=ad_delta- 2*6./delta.^3.*ad_phi_b;
ad_phi_b=0;
%11 tl_delta = ...
%     + fdelta*0.09*tl_p1.*ks./ht ...
%     + fdelta*0.09*p1.*tl_ks./ht ...
%     - fdelta*0.09*p1.*ks./ht.^2.*tl_ht;
ad_p1=ad_p1+ fdelta*0.09.*ks./ht      *ad_delta;
ad_ks=ad_ks+ fdelta*0.09*p1./ht       *ad_delta;
ad_ht=ad_ht- fdelta*0.09*p1.*ks./ht.^2*ad_delta;
ad_delta=0;
%10 tl_p1 = .82*(A./ks).^(.82-1)*( tl_A./ks - A./ks.^2.*tl_ks );
ad_A =ad_A + .82*(A./ks).^(.82-1)./ks     *ad_p1;
ad_ks=ad_ks- .82*(A./ks).^(.82-1)*A./ks.^2*ad_p1;
ad_p1=0;
%9 tl_ht = tl_h - tl_Hm0/2;
ad_h  =ad_h  + 1  *ad_ht;
ad_Hm0=ad_Hm0- 1/2*ad_ht;
ad_ht=0;
%8 tl_Hm0 = tl_Hrms*1.4;
ad_Hrms=ad_Hrms+ 1.4*ad_Hm0;
ad_Hm0=0;
%7 tl_nubar_tb = ...
%     + 2*fw.*uorb.^2/(4*omega).*tl_fw ...
%     + 2*fw.^2.*uorb/(4*omega).*tl_uorb ...
%     - fw.^2.*uorb.^2/(4*omega^2)*tl_omega;
ad_fw   =ad_fw   + 2*fw.*uorb.^2/(4*omega)   *ad_nubar_tb;
ad_uorb =ad_uorb + 2*fw.^2.*uorb/(4*omega)   *ad_nubar_tb;
ad_omega=ad_omega- fw.^2.*uorb.^2/(4*omega^2)*ad_nubar_tb;
ad_nubar_tb=0;
%6 tl_Df = 1/(2*sqrt(pi))*rho*( tl_fw.*uorb.^3 + 3*fw.*uorb.^2.*tl_uorb );
ad_fw  =ad_fw  + 1/(2*sqrt(pi))*rho.*uorb.^3      *ad_Df;
ad_uorb=ad_uorb+ 1/(2*sqrt(pi))*rho*3*fw.*uorb.^2.*ad_Df;
ad_Df=0;
%5 tl_fw = -.52*1.39*(A./z0).^(-.52-1).*( tl_A./z0 - A./z0.^2.*tl_z0 );
ad_A =ad_A - .52*1.39*(A./z0).^(-.52-1)./z0      .*ad_fw;
ad_z0=ad_z0+ .52*1.39*(A./z0).^(-.52-1).*A./z0.^2.*ad_fw;
ad_fw=0;
%4 tl_z0 = tl_ks/33;
ad_ks=ad_ks+1/33*ad_z0;
ad_z0=0;
%3 tl_uorb = tl_A*omega + A*tl_omega;
ad_A    =ad_A    + omega*ad_uorb;
ad_omega=ad_omega+ A    *ad_uorb;
ad_uorb=0;
%2 tl_A = tl_Hrms./(2*sinh(kabs.*h)) ...
%        - Hrms./(2*sinh(kabs.*h)).^2.*( 2*cosh(kabs.*h).*( tl_kabs.*h + kabs.*tl_h ) );
ad_Hrms=ad_Hrms+ 1./(2*sinh(kabs.*h))                              *ad_A;
ad_kabs=ad_kabs- Hrms./(2*sinh(kabs.*h)).^2.*2*cosh(kabs.*h).*h   .*ad_A;
ad_h   =ad_h   - Hrms./(2*sinh(kabs.*h)).^2.*2*cosh(kabs.*h).*kabs.*ad_A;
ad_A=0;
%1 tl_kabs = 1./sqrt(k(:,1).^2+k(:,2).^2).*(k(:,1).*tl_k(:,1)+k(:,2).*tl_k(:,2));
ad_k(:,1)=ad_k(:,1)+ 1./sqrt(k(:,1).^2+k(:,2).^2).*k(:,1).*ad_kabs;
ad_k(:,2)=ad_k(:,2)+ 1./sqrt(k(:,1).^2+k(:,2).^2).*k(:,2).*ad_kabs;
ad_kabs=0;
