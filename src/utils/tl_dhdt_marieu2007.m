function [tl_hpout,tl_qpout]=tl_dhdt_marieu2007(tl_qin,tl_hin,bkgd)

% unpack background struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} '=bkgd.' fld{i} ';']);
end
beta=bkgd.beta;

% Add ghost points
tl_h=zeros(nx,1);  % init
if(ghost==+1)
  tl_q(1:nxin)=tl_qin;
  tl_q(nxin+1)=tl_qin(nxin);
  tl_h(1:nxin)=-tl_hin;
  tl_h(nxin+1)=-tl_hin(nxin);
else
  tl_q(1)=tl_qin(1);
  tl_q(2:(nxin+1))=tl_qin;
  tl_h(1)=-tl_hin(1);
  tl_h(2:(nxin+1))=-tl_hin;
end

% local bed celerity a, eqn (21), centered difference valid pts only
for i=i_nointerp(:)'
  if(i==1)
    tl_a(i) = (tl_q(i+1)-tl_q(i))/(h(i+1)-h(i)) ...
              - (q(i+1)-q(i))/(h(i+1)-h(i))^2*(tl_h(i+1)-tl_h(i));
  elseif(i==nx)
    tl_a(i) = (tl_q(i)-tl_q(i-1))/(h(i)-h(i-1)) ...
              - (q(i)-q(i-1))/(h(i)-h(i-1))^2*(tl_h(i)-tl_h(i-1));
  else
    tl_a(i) = 2*(tl_q(i+1)-tl_q(i-1))/(h(i+1)-h(i-1)) ...
              - 2*(q(i+1)-q(i-1))/(h(i+1)-h(i-1))^2*(tl_h(i+1)-tl_h(i-1));
  end
end

% v2: Replace the potentially-bad points with 3rd order polynomial.  This
% version uses polynomial regression written out in matrix form.  This is
% verified to give the same result as v1 above (polyval/polyfit), but makes
% it easier to write the TL code.  Commented out of NL code after verifying
% that it matches the polyval/polyfit version.  Verification used duck94
% offshore bar migration case (Oct 4, 1994).
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
  tl_pa = Xi*tl_a(ind)';
  tl_a(j)=0;
  for n=0:norder
    tl_a(j) = tl_a(j) + tl_pa(n+1)*j^(norder-n);
  end
end

% interpolate to get q at predictor time step dt/2, eqn 20
tl_xp=-tl_a*dt/2;
tl_qp=zeros(nx,1);
for i=1:nx
  i0=max(find(x<xp(i)));
  if(isempty(i0)) % extrapolate, x<x(1)
    tl_qp(i) = tl_q(1) ...
        + (tl_q(2)-tl_q(1))/dx*(xp(i)-x(1)) ...
        + (q(2)-q(1))/dx*tl_xp(i);
  elseif(i0==nx)  % extrapolate, x>x(nx)
    tl_qp(i) = tl_q(nx) ...
        + (tl_q(nx)-tl_q(nx-1))/dx*(xp(i)-x(nx)) ...
        + (q(nx)-q(nx-1))/dx*tl_xp(i);
  else  % interpolate
    tl_qp(i) = tl_q(i0) ...
        + (tl_q(i0+1)-tl_q(i0))/dx*(xp(i)-x(i0)) ...
        + (q(i0+1)-q(i0))/dx*tl_xp(i);
  end
end

% define candidate functions dh1,dh2,dh3, then apply limiters (MinMod
% function) to get dh, (h' in Marieu et al.), see non-numbered eqns
% following eqn 17
tl_dh123(2:nx,1) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
tl_dh123(1,1) = beta*(tl_h(2)-tl_h(1));
tl_dh123(2:nx-1,2) = .5*(tl_h(3:nx)-tl_h(1:nx-2));
tl_dh123(1,2) = tl_h(2)-tl_h(1);
tl_dh123(nx,2) = tl_h(nx)-tl_h(nx-1);
tl_dh123(1:nx-1,3) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
tl_dh123(nx,3) = beta*(tl_h(nx)-tl_h(nx-1));
tl_dh=zeros(nx,1);
for i=2:nx-1
  if(dh1(i)>0 & dh2(i)>0 & dh3(i)>0)
    [~,ind]=min([dh1(i) dh2(i) dh3(i)]);
    tl_dh(i)=tl_dh123(i,ind);
  elseif(dh1(i)<0 & dh2(i)<0 & dh3(i)<0)
    [~,ind]=max([dh1(i) dh2(i) dh3(i)]);
    tl_dh(i)=tl_dh123(i,ind);
  else
    tl_dh(i)=0;
  end
end
tl_dh(1)=tl_dh(2);
tl_dh(nx)=tl_dh(nx-1);

% h at time step t+dt, eqn 17.  Output is on the staggered grid
tl_hp=zeros(nx-1,1);
tl_hp(1:nx-1) = .5*(tl_h(1:nx-1)+tl_h(2:nx)) ...
    + .125*(tl_dh(1:nx-1)-tl_dh(2:nx)) ...
    - dt/dx*(tl_qp(2:nx)-tl_qp(1:nx-1));

% output is supposed to be on the staggered grid, so interpolate q accordingly
tl_qpout=.5*(tl_qp(1:nx-1)+tl_qp(2:nx));

% revert the result to h=depth instead of h=elevation (see start of code
% where this was flipped)
tl_hpout=-tl_hp;
