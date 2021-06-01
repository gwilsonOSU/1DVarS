function [ad_omega,ad_r,ad_phi] = ad_Uwave_ruessink2012_tcrest(ad_t,omega,r,phi,t)
%
% AD-code for Uwave_ruessink2012_tcrest.m
%
% The value of t should be the point where du/dt = 0 (see
% Uwave_ruessink2012_tcrest.m).  This can be either a maximum or a minimum,
% depending on if you are looking at the crest or the trough.  The code/eqns
% is the same for crest vs. trough.
%

% init ad vars
ad_r=0;
ad_phi=0;
ad_omega=0;

% f,g,h functions in ruessink2012_cresttrough_deriv/.  Here, rename them as
% Ff,Fg,Fh, to avoid conflicts with existing variables names f,g,h
Ff = sin(omega*t) + r*sin(phi)/(1+sqrt(1-r^2));
Fg = r*sin(omega*t + phi);
Fh = 1 - r*cos(omega*t + phi);

% dt/dr
nums1 = Fg./Fh.*( sin(phi)./(1+sqrt(1-r^2)) ...
                  + r^2*sin(phi)/(1+sqrt(1-r^2))^2/sqrt(1-r^2) ) ...
        + Ff/Fh*sin(omega*t+phi) ...
        + Ff*Fg/Fh^2*cos(omega*t+phi);
dens1 = omega*sin(omega*t) + Fg*omega/Fh*cos(omega*t) ...
        + Ff*r*omega/Fh*cos(omega*t+phi) ...
        - Ff*Fg*r*omega/Fh^2*sin(omega*t+phi);
dtdr = -nums1/dens1;

% dt/dphi
nums2 = Fg*r/Fh*cos(phi)/(1+sqrt(1-r^2)) ...
        + Ff*r/Fh*cos(omega*t+phi) ...
        - Ff*Fg*r/Fh^2*sin(omega*t+phi);
dens2 = omega*sin(omega*t) ...
        + Fg*omega/Fh*cos(omega*t) ...
        + Ff*r*omega/Fh*cos(omega*t+phi) ...
        - Ff*Fg*r*omega/Fh^2*sin(omega*t+phi);
dtdphi = -nums2/dens2;

% dt/domega
nums3 = t*sin(omega*t) ...
        + Fg*t/Fh*cos(omega*t) ...
        + Ff/Fh*r*t*cos(omega*t+phi) ...
        - Ff*Fg/Fh^2*r*t*sin(omega*t+phi);
dens3 = omega*sin(omega*t) ...
        + Fg*omega/Fh*cos(omega*t) ...
        + Ff/Fh*r*omega*cos(omega*t+phi) ...
        - Ff*Fg/Fh^2*r*omega*sin(omega*t+phi);
dtdomega = -nums3/dens3;

% tl_t = dtdr*tl_r + dtdphi*tl_phi + dtdomega*tl_omega;
ad_r  =ad_r  + dtdr  *ad_t;
ad_phi=ad_phi+ dtdphi*ad_t;
ad_omega=ad_omega + dtdomega*ad_t;
ad_t=0;

