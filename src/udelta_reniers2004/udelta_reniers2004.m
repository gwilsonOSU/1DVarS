function [udelta,delta_bl,workspc]=udelta_reniers2004(ubar,k,omega,h,Hrms,detady,tau_wind,Dr,fv,ks,d50)%,outvar)
%
% [udelta,delta]=udelta_reniers2004(ubar,k,omega,h,Hrms,detady,tau_wind,Dr,fv,ks,sd50)
%
% Solves for bottom boundary velocity u_delta, see Reniers et al. (2004)
% eqn. (B12).
%
% INPUTS: All must be scalar inputs, unless otherwise stated
%
%   ubar   : depth averaged velocity, 1x2 vector
%   k      : cross-shore wavenumber at spectral peak, rad/m, 1x2 vector
%   omega  : wave frequency at spectral peak, rad/s
%   h      : water depth, m
%   Hrms   : rms wave height, m
%   detady : longshore pressure gradient, m/m units
%   tau_wind: wind stress in N/m2, 1x2 vector
%   Dr     : roller dissipation, W/m2
%   fv     : calibrated, see Reniers et al. (2004) Table 4.  Scalar, of order 0.1
%   d50    : median grain size, m
%

physicalConstants;

% other fixed params
kappa=0.4;  % von karman
fdelta=3;   % default of 3 stated in text (after eqn 16)
nintegral=100;  % number of gridpionts for sigma-integration

% derived params.  See notes here:
% ./reniers2004_solver_notes/reniers2004_formulaTree.sla
kabs=sqrt(k(:,1).^2+k(:,2).^2);
A=Hrms./(2*sinh(kabs.*h));
uorb=A*omega;
z0=ks/33;
fw=1.39*(A./z0).^(-.52);  % eqn (19)
Df=1/(2*sqrt(pi))*rho*fw.*uorb.^3;  % eqn (18)
nubar_tb=fw.^2.*uorb.^2/(4*omega);  % eqn (21)
Hm0=Hrms*1.4;
ht=h-Hm0/2;  % eqn (3)
% delta=fdelta*0.09*(A/ks).^(.82).*ks./ht;  % eqn (15)
p1=(A./ks).^(.82);
delta=fdelta*0.09*p1.*ks./ht;  % eqn (15)
delta_bl = delta.*ht;  % dimensional WBL height
phi_b=6./delta.^2;  % eqn (A9)
tau_wave=Dr.*k/omega;  % eqn (5)
tau_t=tau_wind+tau_wave;  % stated in text
nubar_tflow = 1/6*kappa*ht.*sqrt(g*ht.*abs(detady));  % eqn (11)
abs_tau_wind = sqrt( tau_wind(:,1).^2 + tau_wind(:,2).^2 );
nubar_twind = 1/3*kappa*ht.*sqrt(abs_tau_wind/rho);  % eqn (10)
nubar_twave = fv*Hrms.*(Dr/rho).^(1/3);  % eqn (9)
nubar_t = sqrt( nubar_tflow.^2 ...
                + nubar_twind.^2 ...
                + nubar_twave.^2 );  % eqn (12)
nu_tsurf = 3/2*sqrt( nubar_twind.^2 + nubar_twave.^2 );  % eqn (A5)
sig_s = ( nubar_t - 1/3*nu_tsurf ) ...
        ./( nubar_t - 1/2*nu_tsurf );  % eqn (A4)
phi_s = 1./(sig_s/2-1/3);  % eqn (A3)
nums=phi_s.*nubar_t.*sig_s + phi_b.*nubar_tb.*delta;
dens=phi_s.*nubar_t + phi_b.*nubar_tb;
sig_b = nums./dens;  % eqn (A7)
sig_0 = z0./ht;  % level at which u=0

%--------------------------------------
% Solve for velocity at top of boundary layer, u_delta.  Refer to notes and
% derivations here: ./FgivenUbar_deriv/
%
% Need to do this for each gridpoint separately
%--------------------------------------

% define grid for integration over sigma coordinate
sgridb=linspace(sig_0,delta,nintegral);
sgridm=linspace(delta,1,nintegral);

l1=log(sgridb./sig_0);
l2=log((sig_b-sgridb)./(sig_b-sig_0));
l3=log(sgridm./delta);
l4=log((sig_s-sgridm)./(sig_s-delta));
l1d=log(delta./sig_0);
l2d=log((sig_b-delta)./(sig_b-sig_0));

for i=1:2  % for each direction

% coef functions for bottom boundary layer velocity profile, see
% ./FgivenUbar_deriv/part1_ubot/
Ab(i)=ht./(fv*rho.*phi_s.*nubar_t+rho*phi_b.*nubar_tb);  % eqn (B9)
Bbp(i)=tau_t(i)+Df.*k(i)/omega;
Cbp(i)=-Df.*k(i)./(delta*omega);
alphab=Ab(i).*( Bbp(i)./sig_b.*(l1-l2) - Cbp(i).*l2 );
betab=Ab(i).*( -1/sig_b.*(l1-l2) - l2 );

% coef functions for mid-layer velocity profile, see
% ./FgivenUbar_deriv/part2_umid/.  Note these include terms that
% involve evaluation of alphab,betab at sigma=delta.
Am(i)=ht./(rho*phi_s.*nubar_t);
alpham=alphab(end) ...
       + Am(i).*( tau_t(i)./sig_s.*l3 - tau_t(i)./sig_s.*l4 );
betam=betab(end) ...
      + Am(i).*( -1./sig_s.*l3 + (1./sig_s-1).*l4 );

% define alpha_bar, beta_bar, coefficients for solution ubar.  See 
% ./FgivenUbar_deriv/part3_ubar/
dsb=diff(sgridb(1:2));
dsm=diff(sgridm(1:2));
alphab_bar(i) = trapz(sgridb,alphab);
betab_bar(i)  = trapz(sgridb,betab);
alpham_bar(i) = trapz(sgridm,alpham);
betam_bar(i)  = trapz(sgridm,betam);
alpha_bar(i)  = alphab_bar(i) + alpham_bar(i);
beta_bar(i)   = betab_bar(i) + betam_bar(i);

% solve for F
F(i) = (ubar(i) - alpha_bar(i))./beta_bar(i);

% evaluate eqns (B9)-(B12) at sigma=delta to get u_delta
Bb(i) = tau_t(i) - F(i) + Df.*k(i)/omega;  % eqn (B10)
Cb(i)=F(i)-Df.*k(i)./(delta.*omega);  % eqn (B11)

udelta(i)=Ab(i).*( Bb(i)./sig_b.*l1d ...
            - (Bb(i)./sig_b+Cb(i)).*l2d );  % eqn (B12)

end  % loop on direction i

% % TEST: double-check, should integrate out correctly.  This works.
% ub_fn = @(s)Ab(i)*( Bb(i)/sig_b*log(s/sig_0) ...
%                  - (Bb(i)/sig_b+Cb(i))*log((sig_b-s)/(sig_b-sig_0)) );
% B=tau_t-F;
% C=F;
% um_fn=@(s)ub_fn(delta) + Am(i)*( B/sig_s*log(s/delta) ...
%                              - (B/sig_s+C)*log((sig_s-s)/(sig_s-delta)) );
% ubar_model = integral(ub_fn,sig_0,delta) + integral(um_fn,delta,1);
% if(abs(ubar_model-ubar)>1e-3)
%   disp('ERROR: should never happen!  Dropping to debug mode.')
%   keyboard;
% end

% now save all relevant variables in a struct, so they can be reused in
% TL-AD functions.  Note this works faster if variables are hard-coded.
if(nargout>1)
  workspc=struct;
  workspc.ubar             =ubar             ;
  workspc.k                =k                ;
  workspc.omega            =omega            ;
  workspc.h                =h                ;
  workspc.Hrms             =Hrms             ;
  workspc.detady           =detady           ;
  workspc.Dr               =Dr               ;
  workspc.fv               =fv               ;
  workspc.kabs             =kabs             ;
  workspc.A                =A                ;
  workspc.uorb             =uorb             ;
  workspc.ks               =ks               ;
  workspc.z0               =z0               ;
  workspc.fw               =fw               ;
  workspc.Df               =Df               ;
  workspc.nubar_tb         =nubar_tb         ;
  workspc.Hm0              =Hm0              ;
  workspc.ht               =ht               ;
  workspc.p1               =p1               ;
  workspc.delta            =delta            ;
  workspc.delta_bl         =delta_bl         ;
  workspc.phi_b            =phi_b            ;
  workspc.tau_wind         =tau_wind         ;
  workspc.tau_wind         =tau_wind         ;
  workspc.tau_wave         =tau_wave         ;
  workspc.tau_t            =tau_t            ;
  workspc.nubar_tflow      =nubar_tflow      ;
  workspc.abs_tau_wind     =abs_tau_wind     ;
  workspc.nubar_twind      =nubar_twind      ;
  workspc.nubar_twave      =nubar_twave      ;
  workspc.nubar_t          =nubar_t          ;
  workspc.nu_tsurf         =nu_tsurf         ;
  workspc.sig_s            =sig_s            ;
  workspc.phi_s            =phi_s            ;
  workspc.nums             =nums             ;
  workspc.dens             =dens             ;
  workspc.sig_b            =sig_b            ;
  workspc.sig_0            =sig_0            ;
  workspc.l1d              =l1d              ;
  workspc.l2d              =l2d              ;
  workspc.Ab               =Ab               ;
  workspc.Bbp              =Bbp              ;
  workspc.Cbp              =Cbp              ;
  workspc.Am               =Am               ;
  workspc.alphab_bar       =alphab_bar       ;
  workspc.betab_bar        =betab_bar        ;
  workspc.alpham_bar       =alpham_bar       ;
  workspc.betam_bar        =betam_bar        ;
  workspc.alpha_bar        =alpha_bar        ;
  workspc.beta_bar         =beta_bar         ;
  workspc.F                =F                ;
  workspc.Bb               =Bb               ;
  workspc.Cb               =Cb               ;
  workspc.udelta           =udelta           ;
end

% TEST
% eval(['udelta=' outvar ';']);
