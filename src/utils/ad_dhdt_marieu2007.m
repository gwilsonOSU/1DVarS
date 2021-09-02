function [ad_qin,ad_hin]=ad_dhdt_marieu2007(ad_hpout,ad_qpout,bkgd)

% unpack NL background info
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} '=bkgd.' fld{i} ';']);
end
beta=bkgd.beta;

%-----------------------------------------
% init AD
%-----------------------------------------

ad_qp=zeros(nx,1);
ad_h=zeros(nx,1);
ad_hin=zeros(nxin,1);
ad_qin=zeros(nxin,1);
ad_dh=zeros(nx,1);
ad_dh1=zeros(nx,1);
ad_dh2=zeros(nx,1);
ad_dh3=zeros(nx,1);
ad_q=zeros(nx,1);
ad_xp=zeros(nx,1);
ad_a=zeros(nx,1);
ad_a_orig=zeros(nx,1);

%-----------------------------------------
% begin AD
%-----------------------------------------

%9 revert the result to h=depth instead of h=elevation (see start of code
% where this was flipped)
%9 tl_hpout=-tl_hp;
ad_hp = - ad_hpout;
ad_hpout=0;

%8 output is supposed to be on the staggered grid, so interpolate q accordingly
%8 tl_qpout=.5*(tl_qp(1:nx-1)+tl_qp(2:nx));
ad_qp(1:nx-1)= ad_qp(1:nx-1) + .5*ad_qpout;
ad_qp(2:nx)  = ad_qp(2:nx)+ .5*ad_qpout;
ad_qpout=0;

%7 h at time step t+dt, eqn 17.  Output is on the staggered grid
%7b tl_hp(1:nx-1) = .5*(tl_h(1:nx-1)+tl_h(2:nx)) ...
%     + .125*(tl_dh(1:nx-1)-tl_dh(2:nx)) ...
%     - dt/dx*(tl_qp(2:nx)-tl_qp(1:nx-1));
ad_h(1:nx-1) =ad_h(1:nx-1) + .5   *ad_hp(1:nx-1);
ad_h(2:nx)   =ad_h(2:nx)   + .5   *ad_hp(1:nx-1);
ad_dh(1:nx-1)=ad_dh(1:nx-1)+ .125 *ad_hp(1:nx-1);
ad_dh(2:nx)  =ad_dh(2:nx)  - .125 *ad_hp(1:nx-1);
ad_qp(2:nx)  =ad_qp(2:nx)  - dt/dx*ad_hp(1:nx-1);
ad_qp(1:nx-1)=ad_qp(1:nx-1)+ dt/dx*ad_hp(1:nx-1);
ad_hp(1:nx-1)=0;
%7a tl_hp=nan(nx-1,1);
ad_hp=nan(nx-1,1);

%6 define candidate functions dh1,dh2,dh3, then apply limiters (MinMod
% function) to get dh, (h' in Marieu et al.), see non-numbered eqns
% following eqn 17

%6k tl_dh(nx)=tl_dh(nx-1);
ad_dh(nx-1) = ad_dh(nx-1) + ad_dh(nx);
ad_dh(nx)=0;

%6j tl_dh(1)=tl_dh(2);

% 6i
for i=nx-1:-1:2
  if((dh1(i)>0 & dh2(i)>0 & dh3(i)>0) || (dh1(i)<0 & dh2(i)<0 & dh3(i)<0))
    if(dh(i)==dh1(i))
      %6i1 tl_dh(i)=tl_dh1(i);
      ad_dh1(i)=ad_dh1(i)+ad_dh(i);
      ad_dh(i)=0;
    elseif(dh(i)==dh2(i))
      %6i2 tl_dh(i)=tl_dh2(i);
      ad_dh2(i)=ad_dh2(i)+ad_dh(i);
      ad_dh(i)=0;
    else
      %6i3 tl_dh(i)=tl_dh3(i);
      ad_dh3(i)=ad_dh3(i)+ad_dh(i);
      ad_dh(i)=0;
    end
  else
    %6i4 tl_dh(i)=0;
    ad_dh(i)=0;
  end
end

%6h tl_dh=nan(nx,1);
ad_dh=nan(nx,1);

%6g tl_dh3(nx) = beta*(tl_h(nx)-tl_h(nx-1));
ad_h(nx)  =ad_h(nx)  + beta*ad_dh3(nx);
ad_h(nx-1)=ad_h(nx-1)- beta*ad_dh3(nx);
ad_dh3(nx)=0;

%6f tl_dh3(1:nx-1) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
ad_h(2:nx)  =ad_h(2:nx)  + beta*ad_dh3(1:nx-1);
ad_h(1:nx-1)=ad_h(1:nx-1)- beta*ad_dh3(1:nx-1);
ad_dh3(1:nx-1)=0;

%6e tl_dh2(nx) = tl_h(nx)-tl_h(nx-1);
ad_h(nx)  =ad_h(nx)  + ad_dh2(nx);
ad_h(nx-1)=ad_h(nx-1)- ad_dh2(nx);
ad_dh2(nx)=0;

%6d tl_dh2(1) = tl_h(2)-tl_h(1);
ad_h(2)=ad_h(2)+ ad_dh2(1);
ad_h(1)=ad_h(1)- ad_dh2(1);
ad_dh2(1)=0;

%6c tl_dh2(2:nx-1) = .5*(tl_h(3:nx)-tl_h(1:nx-2));
ad_h(3:nx)  =ad_h(3:nx)  + .5*ad_dh2(2:nx-1);
ad_h(1:nx-2)=ad_h(1:nx-2)- .5*ad_dh2(2:nx-1);
ad_dh2(2:nx-1)=0;

%6b tl_dh1(1) = beta*(tl_h(2)-tl_h(1));
ad_h(2)=ad_h(2)+beta*ad_dh1(1);
ad_h(1)=ad_h(1)-beta*ad_dh1(1);
ad_dh1(1)=0;

%6a tl_dh1(2:nx) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
ad_h(2:nx)  =ad_h(2:nx)  + beta*ad_dh1(2:nx);
ad_h(1:nx-1)=ad_h(1:nx-1)- beta*ad_dh1(2:nx);
ad_dh1(2:nx)=0;

%5 interpolate to get q at predictor time step dt/2, eqn 20
for i=nx:-1:1
  i0=max(find(x<xp(i)));
  ad_q0=0;  % init
  ad_q1=0;  % init
  if(isempty(i0)) % extrapolate, x<x(1)
    x0=x(1)-dx;
    q0=2*q(2)-q(1);
    %5d1 tl_qp(i) = tl_q0 + (tl_q(1)-tl_q0)/dx*(xp(i)-x0) ...
    %     + (q(1)-q0)/dx*tl_xp(i);
    ad_q0   =ad_q0   +               1*ad_qp(i);  
    ad_q(1) =ad_q(1) + 1/dx*(xp(i)-x0)*ad_qp(i);
    ad_q0   =ad_q0   - 1/dx*(xp(i)-x0)*ad_qp(i);
    ad_xp(i)=ad_xp(i)+ (q(1)-q0)/dx   *ad_qp(i);
    ad_qp(i)=0;
    %5c1 tl_q0=2*tl_q(2)-tl_q(1);
    ad_q(2)=ad_q(2)+ 2*ad_q0;
    ad_q(1)=ad_q(1)- 1*ad_q0;
    ad_q0=0;
  elseif(i0==nx)  % extrapolate, x>x(nx)
    x1=x(i0)+dx;
    q1=2*q(i0)-q(i0-1);
    %5d2 tl_qp(i) = tl_q(i0) + (tl_q1-tl_q(i0))/dx*(xp(i)-x(i0)) ...
    %     + (q1-q(i0))/dx*tl_xp(i);
    tl_q(i0)=tl_q(i0)+                  1*ad_qp(i);
    tl_q1   =tl_q1   + 1/dx*(xp(i)-x(i0))*ad_qp(i);
    tl_q(i0)=tl_q(i0)- 1/dx*(xp(i)-x(i0))*ad_qp(i);
    tl_xp(i)=tl_xp(i)+ (q1-q(i0))/dx     *ad_qp(i);
    ad_qp(i)=0;
    %5c2 tl_q1=2*tl_q(i0)-tl_q(i0-1);
    ad_q(i0)  =ad_q(i0)  + 2*ad_q1;
    ad_q(i0-1)=ad_q(i0-1)- 1*ad_q1;
    ad_q1=0;
  else  % interpolate
    %5c3 tl_qp(i) = tl_q(i0) + (tl_q(i0+1)-tl_q(i0))/dx*(xp(i)-x(i0)) ...
    %     + (q(i0+1)-q(i0))/dx*tl_xp(i);
    ad_q(i0)  =ad_q(i0)  +                  1*ad_qp(i);
    ad_q(i0+1)=ad_q(i0+1)+ 1/dx*(xp(i)-x(i0))*ad_qp(i);
    ad_q(i0)  =ad_q(i0)  - 1/dx*(xp(i)-x(i0))*ad_qp(i);
    ad_xp(i)  =ad_xp(i)  + (q(i0+1)-q(i0))/dx*ad_qp(i);
    ad_qp(i)=0;
  end
end
%5b tl_qp=nan(nx,1);
ad_qp=nan(nx,1);
%5a tl_xp=-tl_a*dt/2;
ad_a = ad_a - dt/2*ad_xp;
ad_xp=0;

%4 replace with 3rd order polynomial in cases where h has low slope, or where
% there is a local min/max in h
warning off
for j=i_interp(:)'
  ind=unique(min(nx,max(1,j+[-2 -1 +1 +2])));
  ind=setdiff(ind,j);
  norder=length(ind)-1;
  ad_pa=zeros(norder+1,1);
  X=nan(length(ind));
  for n=0:norder
    X(:,norder-n+1)=ind.^n;
  end
  Xi=pinv(X);  % note Xi is constant, not TL-dependent
  for n=norder:-1:0
    %4c tl_a(j) = tl_a(j) + tl_pa(n+1)*j^(norder-n);
    ad_a(j)   =ad_a(j)   +            1*ad_a(j);
    ad_pa(n+1)=ad_pa(n+1)+ j^(norder-n)*ad_a(j);
    ad_a(j)=0;
  end
  %4b tl_a(j)=0;
  ad_a(j)=0;
  %4a tl_pa = Xi*tl_a_orig(ind)';
  ad_a_orig(ind) = ad_a_orig(ind) + Xi'*ad_pa;
  ad_pa=0;
end
warning on

%3 local bed celerity a, eqn (21), using centered difference
%3d tl_a_orig=tl_a;
ad_a = ad_a + ad_a_orig;
ad_a_orig=0;

%3c tl_a(nx) = (tl_q(nx)-tl_q(nx-1))/(h(nx)-h(nx-1)) ...
%     - (q(nx)-q(nx-1))/(h(nx)-h(nx-1))^2*(tl_h(nx)-tl_h(nx-1));
ad_q(nx)  =ad_q(nx)  + 1/(h(nx)-h(nx-1))                *ad_a(nx);
ad_q(nx-1)=ad_q(nx-1)- 1/(h(nx)-h(nx-1))                *ad_a(nx);
ad_h(nx)  =ad_h(nx)  - (q(nx)-q(nx-1))/(h(nx)-h(nx-1))^2*ad_a(nx);
ad_h(nx-1)=ad_h(nx-1)+ (q(nx)-q(nx-1))/(h(nx)-h(nx-1))^2*ad_a(nx);
ad_a(nx)=0;

%3b tl_a(1) = (tl_q(2)-tl_q(1))/(h(2)-h(1)) ...
%           - (q(2)-q(1))/(h(2)-h(1))^2*(tl_h(2)-tl_h(1));
ad_q(2)=ad_q(2)+ 1/(h(2)-h(1))            *ad_a(1);
ad_q(1)=ad_q(1)- 1/(h(2)-h(1))            *ad_a(1);
ad_h(2)=ad_h(2)- (q(2)-q(1))/(h(2)-h(1))^2*ad_a(1);
ad_h(1)=ad_h(1)+ (q(2)-q(1))/(h(2)-h(1))^2*ad_a(1);
ad_a(1)=0;

%3a tl_a(2:nx-1) = 2*(tl_q(3:nx)-tl_q(1:nx-2))./(h(3:nx)-h(1:nx-2)) ...
%     - 2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2)).^2.*(tl_h(3:nx)-tl_h(1:nx-2));
ad_q(3:nx)  =ad_q(3:nx)  + 2./(h(3:nx)-h(1:nx-2))                       .*ad_a(2:nx-1);
ad_q(1:nx-2)=ad_q(1:nx-2)- 2./(h(3:nx)-h(1:nx-2))                       .*ad_a(2:nx-1);
ad_h(3:nx)  =ad_h(3:nx)  - 2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2)).^2.*ad_a(2:nx-1);
ad_h(1:nx-2)=ad_h(1:nx-2)+ 2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2)).^2.*ad_a(2:nx-1);
ad_a(2:nx-1)=0;

%2 add ghost points
if(ghost==+1)
  %2b1 tl_h=[tl_hin; tl_hin(nxin)];
  ad_hin(1:nxin) = ad_hin(1:nxin) + ad_h(1:nxin);
  ad_hin(nxin) = ad_hin(nxin) + ad_h(nxin+1);
  ad_h=0;
  %2a1 tl_q=[tl_qin; tl_qin(nxin)];
  ad_qin(1:nxin) = ad_qin(1:nxin) + ad_q(1:nxin);
  ad_qin(nxin) = ad_qin(nxin) + ad_q(nxin+1);
  ad_q=0;
else
  %2b2 tl_h=[tl_hin(1); tl_hin];
  ad_hin(1) = ad_hin(1) + ad_h(1);
  ad_hin(1:nxin) = ad_hin(1:nxin) + ad_h(2:nx);
  ad_h=0;
  %2a2 tl_q=[tl_qin(1); tl_qin];
  ad_qin(1) = ad_qin(1) + ad_q(1);
  ad_qin(1:nxin) = ad_qin(1:nxin) + ad_q(2:nx);
  ad_q=0;
end

% reverse sign convention for calculations
%1 tl_h=-tl_hin;
ad_hin = ad_hin - ad_h;
ad_h=0;
