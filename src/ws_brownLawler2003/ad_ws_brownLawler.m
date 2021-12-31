function ad_d50=ad_ws_brownLawler(ad_ws,d50);
%
% AD-code for ws_brownLawler.m
%

physicalConstants;

% recalculate NL
mu=nu*rho;  % dynamic viscosity of water, Pa*s
dstar=d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);  % eqn (24)
p1=18./dstar.^2;
p2=0.898*(0.936*dstar+1)./(dstar+1);
p3=0.317./dstar;
p4=p1.^p2;
p5=p3.^0.449;
p6=p4 + p5;
ustar=p6.^(-1.114);  % eqn (33)
ws=ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);  % eqn (24)

% init AD
ad_ustar=0;
ad_dstar=0;
ad_d50=0;
ad_p1=0;
ad_p2=0;
ad_p3=0;
ad_p4=0;
ad_p5=0;
ad_p6=0;

%9 tl_ws = tl_ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);
ad_ustar=ad_ustar+ 1./(rho^2./(g*mu*(rhop-rho))).^(1/3)*ad_ws;
ad_ws=0;
%8 tl_ustar = -1.114*p6.^(-1.114-1).*tl_p6;
ad_p6=ad_p6- 1.114*p6.^(-1.114-1).*ad_ustar
ad_ustar=0;
%7 tl_p6 = tl_p4 + tl_p5;
%         - 0.898*(0.936*dstar+1)./(dstar+1).^2.*tl_dstar;
ad_p4   =ad_p4   + 1                                   *ad_p6;
ad_p5   =ad_p5;  + 1                                   *ad_p6;
ad_dstar=ad_dstar- 0.898*(0.936*dstar+1)./(dstar+1).^2.*ad_p6;
ad_p6=0;
%6 tl_p5 = 0.449*p3.^(0.449-1).*tl_p3;
ad_p3=ad_p3+ 0.449*p3.^(0.449-1).*ad_p5
ad_p5=0;
%5 tl_p4 = p2.*p1.^(p2-1).*tl_p1 + p1.^p2.*log(p1).*tl_p2;
ad_p1=ad_p1+ p2.*p1.^(p2-1) .*ad_p4;
ad_p2=ad_p2+ p1.^p2.*log(p1).*ad_p4;
ad_p4=0;
%4 tl_p3 = -0.317./dstar.^2.*tl_dstar;
ad_dstar=ad_dstar- 0.317./dstar.^2.*ad_p3;
ad_p3=0;
%3 tl_p2 = 0.898*(0.936*tl_dstar)./(dstar+1) ...
%         - 0.898*(0.936*dstar+1)./(dstar+1).^2.*tl_dstar;
ad_dstar=ad_dstar+ 0.898*0.936./(dstar+1)             .*ad_p2;
ad_dstar=ad_dstar- 0.898*(0.936*dstar+1)./(dstar+1).^2.*ad_p2;
ad_p2=0;
%2 tl_p1 = -2*18./dstar.^3.*tl_dstar;
ad_dstar=ad_dstar- 2*18./dstar.^3.*ad_p1;
ad_p1=0;
%1 tl_dstar = tl_d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);
ad_d50=ad_d50+(g*rho*(rhop-rho)/mu^2).^(1/3).*ad_dstar;
ad_dstar=0;
