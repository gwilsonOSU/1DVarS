function dqdx=ddx_upwind(x,q,h)
%
% dqdx=ddx_upwind(x,q)
%
% assumes x is monotonically increasing, and no-flux (q=0) through
% boundaries
%
% x and q should both be Nx1 vectors, where N is the number of gridpoints
%

nx=length(x);

% Local bedform celerity a, using centered difference.  This is based on
% Marieu et al. (2007).  They point out that
%
% dh/dt = dq/dx = (dq/dh)*(dh/dx)
%
% Which can be viewed as a wave equation for h, with c = -dq/dh
a(2:nx-1)=-2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2));
a(1)=-(q(2)-q(1))/(h(2)-h(1));
a(nx)=-(q(nx)-q(nx-1))/(h(nx)-h(nx-1));

ip=find(a>=0);
im=find(a< 0);

% add ghost point zeros for no-flux BCs, and adjust indeces accordingly
q=[q(1); q(:); 0];
x=[x(1)-diff(x(1:2)); x(:); x(end)+diff(x(end-1:end))];
ip=ip+1;
im=im+1;

% calculate upwind finite differences
dqdx=nan(nx+2,1);
dqdx(ip)=(q(ip)-q(ip-1))./(x(ip)-x(ip-1));
dqdx(im)=(q(im+1)-q(im))./(x(im+1)-x(im));

% take out ghost points
dqdx=dqdx(2:end-1);
