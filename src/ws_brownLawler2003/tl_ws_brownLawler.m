function tl_ws=tl_ws_brownLawler(tl_d50,d50);
%
% TL-code for ws_brownLawler.m
%

physicalConstants;

mu=nu*rho;  % dynamic viscosity of water, Pa*s
dstar=d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);  % eqn (24)
tl_dstar = tl_d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);
p1=18./dstar.^2;
tl_p1 = -2*18./dstar.^3.*tl_dstar;
p2=0.898*(0.936*dstar+1)./(dstar+1);
tl_p2 = 0.898*(0.936*tl_dstar)./(dstar+1) ...
        - 0.898*(0.936*dstar+1)./(dstar+1).^2.*tl_dstar;
p3=0.317./dstar;
tl_p3 = -0.317./dstar.^2.*tl_dstar;
p4=p1.^p2;
tl_p4 = p2.*p1.^(p2-1).*tl_p1 + p1.^p2.*log(p1).*tl_p2;
p5=p3.^0.449;
tl_p5 = 0.449*p3.^(0.449-1).*tl_p3;
p6=p4 + p5;
tl_p6 = tl_p4 + tl_p5;
ustar=p6.^(-1.114);  % eqn (33)
tl_ustar = -1.114*p6.^(-1.114-1).*tl_p6;

ws=ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);  % eqn (24)
tl_ws = tl_ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);
