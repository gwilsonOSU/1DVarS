function tl_dqdx_out=tl_ddx_centered(tl_q_in,x,q)

nx=length(x);
nxp=nx+2;

% add ghost point zeros for no-flux BCs
tl_q=zeros(nxp,1);  % init
tl_q(1) = tl_q_in(1);
tl_q(2:nx+1) = tl_q_in(1:nx);
tl_q(nx+2) = tl_q_in(nx);
x=[x(1)-diff(x(1:2)); x(:); x(nx)+diff(x(nx-1:nx))];

% centered finite difference
tl_dqdx=zeros(nxp,1);  % init
tl_dqdx(2:nxp-1) = (tl_q(3:nxp)-tl_q(1:nxp-2))./(x(3:nxp)-x(1:nxp-2));

% take out ghost points
tl_dqdx_out=tl_dqdx(2:nxp-1);
