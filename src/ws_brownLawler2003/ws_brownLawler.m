function ws=ws_brownLawler(d50);
%
% ws=ws_brownLawler(d50);
%
% settling velocity for quartz sand, from Brown & Lawler
%

physicalConstants;

mu=nu*rho;  % dynamic viscosity of water, Pa*s
dstar=d50.*(g*rho*(rhop-rho)/mu^2).^(1/3);  % eqn (24)

% % v1: one line equation
% ustar0=( (18./dstar.^2).^(0.898*(0.936*dstar+1)./(dstar+1)) ...
%         + (0.317./dstar).^0.449 ).^(-1.114);  % eqn (33)

% v2: split into multiple lines, this makes TL-AD coding easier.
% Double-checked, this gives the same answer as v1 above
p1=18./dstar.^2;
p2=0.898*(0.936*dstar+1)./(dstar+1);
p3=0.317./dstar;
p4=p1.^p2;
p5=p3.^0.449;
p6=p4 + p5;
ustar=p6.^(-1.114);  % eqn (33)

ws=ustar./(rho^2./(g*mu*(rhop-rho))).^(1/3);  % eqn (24)
