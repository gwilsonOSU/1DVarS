%
% test of derived equations for TL model of crest and trough timing.  See
% ruessink2012_cresttrough_deriv/
%
clear

% input params, somewhat arbitrary
h=2;
Hmo=2;
omega=2*pi/10;
k=fzero(@(k)omega^2-9.8*k*tanh(k*h),.1);
t=linspace(0,2*pi/omega,1000);

% bkgd NL model
[u,wksp]=Uwave_ruessink2012(omega*t,Hmo,k,omega,h);
r  =wksp.r  ;
phi=wksp.phi;

% solve for crest and trough times
dudt=@(t)cos(omega*t) ...
     - (sin(omega*t)+r*sin(phi)/(1+sqrt(1-r^2))) ...
     /(1-r*cos(omega*t+phi)) ...
     *r*sin(omega*t+phi);
[~,icu_guess]=max(u);
[~,itu_guess]=min(u);
Tc=fzero(dudt,t(icu_guess));
Tt=fzero(dudt,t(itu_guess));

% perturb and re-solve the NL model
tl_frac=0.001;
myrnd=@()2*(1-rand(1));
tl_omega=omega*tl_frac*myrnd();
tl_r  =r  *tl_frac*myrnd();
tl_phi=phi*tl_frac*myrnd();
dudt2=@(t)cos((omega+tl_omega)*t) ...
     - (sin((omega+tl_omega)*t)+(r+tl_r)*sin((phi+tl_phi))/(1+sqrt(1-(r+tl_r)^2))) ...
      /(1-(r+tl_r)*cos((omega+tl_omega)*t+(phi+tl_phi))) ...
     *(r+tl_r)*sin((omega+tl_omega)*t+(phi+tl_phi));
Tc2=fzero(dudt2,Tc);
Tt2=fzero(dudt2,Tt);

% TL-model
tl_Tc = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r,tl_phi,omega,r,phi,Tc);
tl_Tt = tl_Uwave_ruessink2012_tcrest(tl_omega,tl_r,tl_phi,omega,r,phi,Tt);

% compare
disp(['tl_Tc: true = ' num2str(Tc2-Tc)])
disp(['       TL   = ' num2str(tl_Tc)])
disp(['tl_Tt: true = ' num2str(Tt2-Tt)])
disp(['       TL   = ' num2str(tl_Tt)])
