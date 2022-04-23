function dqdx=ddx_centered(x,q)
%
% dqdx=ddx_centered(x,q)
%
% x and q should both be Nx1 vectors, where N is the number of gridpoints
%

nx=length(x);

% add ghost point zeros for no-flux BCs
q=[q(1); q(:); q(end)];
x=[x(1)-diff(x(1:2)); x(:); x(end)+diff(x(end-1:end))];

% centered finite difference
dqdx=zeros(nx+2,1);  % init
dqdx(2:end-1)=(q(3:end)-q(1:end-2))./(x(3:end)-x(1:end-2));

% take out ghost points
dqdx=dqdx(2:end-1);
