function tl_dqdx=tl_ddx_upwind(tl_q,x,q,h)

nx=length(x);

a(2:nx-1)=-2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2));
a(1)=-(q(2)-q(1))/(h(2)-h(1));
a(nx)=-(q(nx)-q(nx-1))/(h(nx)-h(nx-1));

% % oops, I don't actually need the TL of the bed celerity
% tl_a(2:nx-1)=-2*(tl_q(3:nx)-tl_q(1:nx-2))./(h(3:nx)-h(1:nx-2)) ...
%     + 2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2)).^2.*(tl_h(3:nx)-tl_h(1:nx-2));
% tl_a(1)=-(tl_q(2)-tl_q(1))/(h(2)-h(1)) ...
%         + (q(2)-q(1))/(h(2)-h(1))^2*(tl_h(2)-tl_h(1));
% tl_a(nx)=-(tl_q(nx)-tl_q(nx-1))/(h(nx)-h(nx-1)) ...
%          + (q(nx)-q(nx-1))/(h(nx)-h(nx-1))^2*(tl_h(nx)-tl_h(nx-1));

ip=find(a>=0);
im=find(a< 0);

% add ghost point zeros for no-flux BCs, and adjust indeces accordingly
qp=[q(1); q(:); 0];
tl_qp=[tl_q(1); tl_q(:); 0];
xp=[x(1)-diff(x(1:2)); x(:); x(end)+diff(x(end-1:end))];
ip=ip+1;
im=im+1;

% calculate upwind finite differences
tl_dqpdx=nan(nx+2,1);
tl_dqpdx(ip)=(tl_qp(ip)-tl_qp(ip-1))./(xp(ip)-xp(ip-1));
tl_dqpdx(im)=(tl_qp(im+1)-tl_qp(im))./(xp(im+1)-xp(im));

% take out ghost points
tl_dqdx=tl_dqpdx(2:end-1);
