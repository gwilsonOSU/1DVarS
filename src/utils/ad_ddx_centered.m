function ad_q_in=ad_ddx_centered(ad_dqdx_out,x,q)

nx=length(x);
x=[x(1)-diff(x(1:2)); x(:); x(nx)+diff(x(nx-1:nx))];
nxp=nx+2;

% init AD
ad_dqdx=zeros(nxp,1);
ad_q=zeros(nxp,1);
ad_q_in=zeros(nx,1);

%4 tl_dqdx_out=tl_dqdx(2:nxp-1);
ad_dqdx(2:nxp-1) = ad_dqdx(2:nxp-1) + ad_dqdx_out;
ad_dqdx_out=0;
%3 tl_dqdx(2:nxp-1) = (tl_q(3:nxp)-tl_q(1:nxp-2))./(x(3:nxp)-x(1:nxp-2));
ad_q(3:nxp)  =ad_q(3:nxp)  + 1./(x(3:nxp)-x(1:nxp-2)).*ad_dqdx(2:nxp-1);
ad_q(1:nxp-2)=ad_q(1:nxp-2)- 1./(x(3:nxp)-x(1:nxp-2)).*ad_dqdx(2:nxp-1);
ad_dqdx(2:nxp-1)=0;
%2 tl_dqdx=zeros(nxp,1);  % init
ad_dqdx=zeros(nxp,1);  % init

%1c tl_q(nx+2) = tl_q_in(nx);
ad_q_in(nx) = ad_q_in(nx) + ad_q(nx+2);
%1b tl_q(2:nx+1) = tl_q_in(1:nx);
ad_q_in(1:nx) = ad_q_in(1:nx) + ad_q(2:nx+1);
%1a tl_q(1) = tl_q_in(1);
ad_q_in(1) = ad_q_in(1) + ad_q(1);
ad_q=0;
