function tl_ws=tl_ws_brownLawler(tl_d50,d50);
%
% TL-code for ws_brownLawler.m
%
% Testing with small perturbations suggests this may not be exact, but I
% can't find any typos.  It may just be the equation is very nonlinear

physicalConstants;

% recalculate NL variables
mu=nu*rho;  % dynamic viscosity of water, Pa*s
dstar=d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);  % eqn (24)
p1=18./dstar.^2;
p2=0.898*(0.936*dstar+1)./(dstar+1);
p3=0.317./dstar;
p4=p1.^p2;
p5=p3.^0.449;
p6=p4 + p5;
% ustar=p6.^(-1.114);  % eqn (33)
% ws=ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);  % eqn (24)

tl_dstar=tl_d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);  % ok
tl_p1=-2*18./dstar.^3.*tl_dstar; % ok
tl_p2=0.898*( 0.936./(dstar+1).*tl_dstar ...
              - (0.936*dstar+1)./(dstar+1).^2.*tl_dstar ); % ok
tl_p3=0.317./dstar.^2.*tl_dstar;  % ok
tl_p4 = p2.*p1.^(p2-1).*tl_p1 ...
        + p1.^p2.*log(p1).*tl_p2;  % ok, double-checked numerically
tl_p5=0.449*p3.^(0.449-1).*tl_p3;  % ok
tl_p6=tl_p4 + tl_p5; % ok
tl_ustar=-1.114*p6.^(-1.114-1).*tl_p6;  % ok
tl_ws=tl_ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);  % ok
