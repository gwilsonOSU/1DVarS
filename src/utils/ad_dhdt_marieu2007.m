function [ad_qin,ad_hin]=ad_dhdt_marieu2007(ad_hpout,ad_qpout,bkgd)

% unpack background struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} '=bkgd.' fld{i} ';']);
end
beta=bkgd.beta;

%------------------------------------
% init AD
%------------------------------------

ad_qin=zeros(nxin,1);
ad_hin=zeros(nxin,1);
ad_hp=zeros(nxin,1);

ad_dh=zeros(nx,1);
ad_qp=zeros(nx,1);
ad_xp=zeros(nx,1);
ad_q=zeros(nx,1);
ad_h=zeros(nx,1);
ad_dh123=zeros(nx,3);
ad_a=zeros(nx,1);

%------------------------------------
% begin AD
%------------------------------------

% revert the result to h=depth instead of h=elevation
%21 tl_hpout=-tl_hp;
ad_hp = ad_hp - ad_hpout;
ad_hpout=0;

% switch back to staggered grid
%20 tl_qpout=.5*(tl_qp(1:nx-1)+tl_qp(2:nx));
ad_qp(1:nx-1)=ad_qp(1:nx-1)+ .5*ad_qpout;
ad_qp(2:nx)  =ad_qp(2:nx)  + .5*ad_qpout;
ad_qpout=0;

% h at time step t+dt, eqn 17.  Output is on the staggered grid
%19 tl_hp(1:nx-1) = .5*(tl_h(1:nx-1)+tl_h(2:nx)) ...
%     + .125*(tl_dh(1:nx-1)-tl_dh(2:nx)) ...
%     - dt/dx*(tl_qp(2:nx)-tl_qp(1:nx-1));
ad_h(1:nx-1) =ad_h(1:nx-1) + .5   *ad_hp(1:nx-1);
ad_h(2:nx)   =ad_h(2:nx)   + .5   *ad_hp(1:nx-1);
ad_dh(1:nx-1)=ad_dh(1:nx-1)+ .125 *ad_hp(1:nx-1);
ad_dh(2:nx)  =ad_dh(2:nx)  - .125 *ad_hp(1:nx-1);
ad_qp(2:nx)  =ad_qp(2:nx)  - dt/dx*ad_hp(1:nx-1);
ad_qp(1:nx-1)=ad_qp(1:nx-1)+ dt/dx*ad_hp(1:nx-1);
ad_hp(1:nx-1)=0;

%18 tl_hp=zeros(nx-1,1);
ad_hp=zeros(nx-1,1);

%17 tl_dh(nx)=tl_dh(nx-1);
ad_dh(nx-1) = ad_dh(nx-1) + ad_dh(nx);
ad_dh(nx)=0;

%16 tl_dh(1)=tl_dh(2);
ad_dh(2) = ad_dh(2) + ad_dh(1);
ad_dh(1)=0;

% define candidate functions dh1,dh2,dh3, then apply limiters (MinMod
% function) to get dh, (h' in Marieu et al.), see non-numbered eqns
% following eqn 17
for i=2:nx-1
  if(dh1(i)>0 & dh2(i)>0 & dh3(i)>0)
    [~,ind]=min([dh1(i) dh2(i) dh3(i)]);
    %15a tl_dh(i)=tl_dh123(i,ind);
    ad_dh123(i,ind) = ad_dh123(i,ind) + ad_dh(i);
    ad_dh(i)=0;
  elseif(dh1(i)<0 & dh2(i)<0 & dh3(i)<0)
    [~,ind]=max([dh1(i) dh2(i) dh3(i)]);
    %15b tl_dh(i)=tl_dh123(i,ind);
    ad_dh123(i,ind) = ad_dh123(i,ind) + ad_dh(i);
    ad_dh(i)=0;
  else
    %15c tl_dh(i)=0;
    ad_dh(i)=0;
  end
end
%14 tl_dh=zeros(nx,1);
ad_dh=zeros(nx,1);

%13 tl_dh123(nx,3) = beta*(tl_h(nx)-tl_h(nx-1));
ad_h(nx)  =ad_h(nx)  + beta*ad_dh123(nx,3);
ad_h(nx-1)=ad_h(nx-1)- beta*ad_dh123(nx,3);
ad_dh123(nx,3)=0;

%12 tl_dh123(1:nx-1,3) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
ad_h(2:nx)  =ad_h(2:nx)  + beta*ad_dh123(1:nx-1,3);
ad_h(1:nx-1)=ad_h(1:nx-1)- beta*ad_dh123(1:nx-1,3);
ad_dh123(1:nx-1,3)=0;

%11 tl_dh123(nx,2) = tl_h(nx)-tl_h(nx-1);
ad_h(nx)  =ad_h(nx)  + ad_dh123(nx,2);
ad_h(nx-1)=ad_h(nx-1)- ad_dh123(nx,2);
ad_dh123(nx,2)=0;

%10 tl_dh123(1,2) = tl_h(2)-tl_h(1);
ad_h(2)=ad_h(2)+ ad_dh123(1,2);
ad_h(1)=ad_h(1)- ad_dh123(1,2);
ad_dh123(1,2)=0;

%9 tl_dh123(2:nx-1,2) = .5*(tl_h(3:nx)-tl_h(1:nx-2));
ad_h(3:nx)  =ad_h(3:nx)  + .5*ad_dh123(2:nx-1,2);
ad_h(1:nx-2)=ad_h(1:nx-2)- .5*ad_dh123(2:nx-1,2);
ad_dh123(2:nx-1,2)=0;

%8 tl_dh123(1,1) = beta*(tl_h(2)-tl_h(1));
ad_h(2)=ad_h(2)+ beta*ad_dh123(1,1);
ad_h(1)=ad_h(1)- beta*ad_dh123(1,1);
ad_dh123(1,1)=0;

%7 tl_dh123(2:nx,1) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
ad_h(2:nx)  =ad_h(2:nx)  + beta*ad_dh123(2:nx,1);
ad_h(1:nx-1)=ad_h(1:nx-1)- beta*ad_dh123(2:nx,1);
ad_dh123(2:nx,1)=0;

% interpolate to get q at predictor time step dt/2, eqn 20
for i=1:nx
  i0=max(find(x<xp(i)));
  if(isempty(i0)) % extrapolate, x<x(1)
    % tl_qp(i) = tl_q(1) ...
    %     + (tl_q(2)-tl_q(1))/dx*(xp(i)-x(1)) ...
    %     + (q(2)-q(1))/dx*tl_xp(i);
    ad_q(1) =ad_q(1) + 1                *ad_qp(i);
    ad_q(2) =ad_q(2) + 1/dx*(xp(i)-x(1))*ad_qp(i);
    ad_q(1) =ad_q(1) - 1/dx*(xp(i)-x(1))*ad_qp(i);
    ad_xp(i)=ad_xp(i)+ (q(2)-q(1))/dx   *ad_qp(i);
    ad_qp(i)=0;
  elseif(i0==nx)  % extrapolate, x>x(nx)
    % tl_qp(i) = tl_q(nx) ...
    %     + (tl_q(nx)-tl_q(nx-1))/dx*(xp(i)-x(nx)) ...
    %     + (q(nx)-q(nx-1))/dx*tl_xp(i);
    ad_q(nx)  =ad_q(nx)  + 1                 *ad_qp(i);
    ad_q(nx)  =ad_q(nx)  + 1/dx*(xp(i)-x(nx))*ad_qp(i);
    ad_q(nx-1)=ad_q(nx-1)- 1/dx*(xp(i)-x(nx))*ad_qp(i);
    ad_xp(i)  =ad_xp(i)  + (q(nx)-q(nx-1))/dx*ad_qp(i);
    ad_qp(i)=0;
  else  % interpolate
    % tl_qp(i) = tl_q(i0) ...
    %     + (tl_q(i0+1)-tl_q(i0))/dx*(xp(i)-x(i0)) ...
    %     + (q(i0+1)-q(i0))/dx*tl_xp(i);
    ad_q(i0)  =ad_q(i0)  + 1                 *ad_qp(i);
    ad_q(i0+1)=ad_q(i0+1)+ 1/dx*(xp(i)-x(i0))*ad_qp(i);
    ad_q(i0)  =ad_q(i0)  - 1/dx*(xp(i)-x(i0))*ad_qp(i);
    ad_xp(i)  =ad_xp(i)  + (q(i0+1)-q(i0))/dx*ad_qp(i);
    ad_qp(i)=0;
  end
end
% tl_qp=zeros(nx,1);
ad_qp=zeros(nx,1);
% tl_xp=-tl_a*dt/2;
ad_a = ad_a - dt/2*ad_xp;
ad_xp=0;

%5 tl_qp=zeros(nx,1);
ad_qp=zeros(nx,1);

%4 tl_xp=-tl_a*dt/2;
ad_a=ad_a-dt/2*ad_xp;
ad_xp=0;

%3 Replace the potentially-bad points with 3rd order polynomial
for j=i_interp(:)'
  ind=setdiff(1:nx,[j i_interp(:)']);
  thisdx=x(j)-x(ind);
  [~,isort]=sort(abs(thisdx),'ascend');
  ind=ind(isort(1:4));  % four nearest valid points
  ind=setdiff(ind,j);
  norder=length(ind)-1;
  X=zeros(length(ind));
  for n=0:norder
    X(:,norder-n+1)=ind.^n;
  end
  Xi=pinv(X);  % note Xi is constant, not TL-dependent
  ad_pa=zeros(norder+1,1);  % init
  for n=norder:-1:0
    %3a3 tl_a(j) = tl_a(j) + tl_pa(n+1)*j^(norder-n);
    ad_a(j)   =ad_a(j)   + 1           *ad_a(j);
    ad_pa(n+1)=ad_pa(n+1)+ j^(norder-n)*ad_a(j);
    ad_a(j)=0;
  end
  %3a2 tl_a(j)=0;
  ad_a(j)=0;
  %3a1 tl_pa = Xi*tl_a(ind)';
  ad_a(ind) = ad_a(ind) + Xi'*ad_pa;
  ad_pa=0;
end

%2 local bed celerity a, eqn (21), centered difference valid pts only
for i=i_nointerp(:)'
  if(i==1)
    %2a1 tl_a(i) = (tl_q(i+1)-tl_q(i))/(h(i+1)-h(i)) ...
    %           - (q(i+1)-q(i))/(h(i+1)-h(i))^2*(tl_h(i+1)-tl_h(i));
    ad_q(i+1)=ad_q(i+1)+ 1/(h(i+1)-h(i))              *ad_a(i);
    ad_q(i)  =ad_q(i)  - 1/(h(i+1)-h(i))              *ad_a(i);
    ad_h(i+1)=ad_h(i+1)- (q(i+1)-q(i))/(h(i+1)-h(i))^2*ad_a(i);
    ad_h(i)  =ad_h(i)  + (q(i+1)-q(i))/(h(i+1)-h(i))^2*ad_a(i);
    ad_a(i)=0;
  elseif(i==nx)
    %2b1 tl_a(i) = (tl_q(i)-tl_q(i-1))/(h(i)-h(i-1)) ...
    %           - (q(i)-q(i-1))/(h(i)-h(i-1))^2*(tl_h(i)-tl_h(i-1));
    ad_q(i)  =ad_q(i)  + 1/(h(i)-h(i-1))              *ad_a(i);
    ad_q(i-1)=ad_q(i-1)- 1/(h(i)-h(i-1))              *ad_a(i);
    ad_h(i)  =ad_h(i)  - (q(i)-q(i-1))/(h(i)-h(i-1))^2*ad_a(i);
    ad_h(i-1)=ad_h(i-1)+ (q(i)-q(i-1))/(h(i)-h(i-1))^2*ad_a(i);
    ad_a(i)=0;
  else
    %2c1 tl_a(i) = 2*(tl_q(i+1)-tl_q(i-1))/(h(i+1)-h(i-1)) ...
    %           - 2*(q(i+1)-q(i-1))/(h(i+1)-h(i-1))^2*(tl_h(i+1)-tl_h(i-1));
    ad_q(i+1)=ad_q(i+1)+ 2/(h(i+1)-h(i-1))                  *ad_a(i);
    ad_q(i-1)=ad_q(i-1)- 1/(h(i+1)-h(i-1))                  *ad_a(i);
    ad_h(i+1)=ad_h(i+1)- 2*(q(i+1)-q(i-1))/(h(i+1)-h(i-1))^2*ad_a(i);
    ad_h(i-1)=ad_h(i-1)+ 2*(q(i+1)-q(i-1))/(h(i+1)-h(i-1))^2*ad_a(i);
    ad_a(i)=0;
  end
end

%1 Add ghost points
if(ghost==+1)
  %1a4 tl_h(nxin+1)=-tl_hin(nxin);
  ad_hin(nxin)=ad_hin(nxin)-ad_h(nxin+1);
  ad_h(nxin+1)=0;
  %1a3 tl_h(1:nxin)=-tl_hin;
  ad_hin = ad_hin - ad_h(1:nxin);
  ad_h(1:nxin)=0;
  %1a2 tl_q(nxin+1)=tl_qin(nxin);
  ad_qin(nxin) = ad_qin(nxin) + ad_q(nxin+1);
  ad_q(nxin+1)=0;
  %1a1 tl_q(1:nxin)=tl_qin;
  ad_qin=ad_qin+ad_q(1:nxin);
  ad_q(1:nxin)=0;
else
  %1b4 tl_h(2:(nxin+1))=-tl_hin;
  ad_hin=ad_hin-ad_h(2:(nxin+1));
  ad_h(2:(nxin+1))=0;
  %1b3 tl_h(1)=-tl_hin(1);
  ad_hin(1)=ad_hin(1)-ad_h(1);
  ad_h(1)=0;
  %1b2 tl_q(2:(nxin+1))=tl_qin;
  ad_qin=ad_qin+ad_q(2:(nxin+1));
  ad_q(2:(nxin+1))=0;
  %1b1 tl_q(1)=tl_qin(1);
  ad_qin(1)=ad_qin(1)+ad_q(1);
  ad_q(1)=0;
end
