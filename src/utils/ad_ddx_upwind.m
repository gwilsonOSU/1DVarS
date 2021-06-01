function ad_q=ad_ddx_upwind(ad_dqdx,x,q,h)

nx=length(x);

a(2:nx-1)=-2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2));
a(1)=-(q(2)-q(1))/(h(2)-h(1));
a(nx)=-(q(nx)-q(nx-1))/(h(nx)-h(nx-1));

ip=find(a>=0);
im=find(a< 0);

% add ghost point zeros for no-flux BCs, and adjust indeces accordingly
xp=[x(1)-diff(x(1:2)); x(:); x(end)+diff(x(end-1:end))];
ip=ip+1;
im=im+1;

% init adjoint
ad_q=zeros(nx,1);
ad_qp=zeros(nx+2,1);
ad_dqpdx=zeros(nx+2,1);

%4 tl_dqdx=tl_dqpdx(2:end-1);
%         =[zeros(nx,1) eye(nx) zeros(nx,1)]*tl_dqpdx;
ad_dqpdx=ad_dqpdx+[zeros(nx,1) eye(nx) zeros(nx,1)]'*ad_dqdx;
ad_dqdx=0*ad_dqdx;

%3 tl_dqpdx(ip)=(tl_qp(ip)-tl_qp(ip-1))./(xp(ip)-xp(ip-1));
ad_qp(ip)  =ad_qp(ip)  + 1./(xp(ip)-xp(ip-1)).*ad_dqpdx(ip); 
ad_qp(ip-1)=ad_qp(ip-1)- 1./(xp(ip)-xp(ip-1)).*ad_dqpdx(ip);
ad_dqpdx(ip)=0;

%2 tl_dqpdx(im)=(tl_qp(im+1)-tl_qp(im))./(xp(im+1)-xp(im));
ad_qp(im+1)=ad_qp(im+1)+ 1./(xp(im+1)-xp(im)).*ad_dqpdx(im);
ad_qp(im)  =ad_qp(im)  - 1./(xp(im+1)-xp(im)).*ad_dqpdx(im);
ad_dqpdx(im)=0;

%1 tl_qp=[tl_q(1); tl_q; 0];
%       =[[1 zeros(1,nx-1)]; eye(nx); zeros(1,nx)]*tl_q;
ad_q = ad_q + [[1 zeros(1,nx-1)]; eye(nx); zeros(1,nx)]'*ad_qp;
ad_qp=0*ad_qp;
