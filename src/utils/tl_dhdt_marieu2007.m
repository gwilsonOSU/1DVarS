function [tl_hpout,tl_qpout]=tl_dhdt_marieu2007(tl_qin,tl_hin,bkgd)

% unpack NL background info
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} '=bkgd.' fld{i} ';']);
end

% reverse sign convention for calculations
tl_h=-tl_hin;

% add ghost points
nx=length(bkgd.qin);
if(ghost==+1)
  tl_q=[tl_qin; tl_qin(nx)];
  tl_h=[tl_hin; tl_hin(nx)];
else
  tl_q=[tl_qin(1); tl_qin];
  tl_h=[tl_hin(1); tl_hin];
end

% local bed celerity a, eqn (21), using centered difference
tl_a(2:nx-1) = 2*(tl_q(3:nx)-tl_q(1:nx-2))./(h(3:nx)-h(1:nx-2)) ...
    - 2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2)).^2.*(tl_h(3:nx)-tl_h(1:nx-2));
tl_a(1) = (tl_q(2)-tl_q(1))/(h(2)-h(1)) ...
          - (q(2)-q(1))/(h(2)-h(1))^2*(tl_h(2)-tl_h(1));
tl_a(nx) = (tl_q(nx)-tl_q(nx-1))/(h(nx)-h(nx-1)) ...
    - (q(nx)-q(nx-1))/(h(nx)-h(nx-1))^2*(tl_h(nx)-tl_h(nx-1));

% replace with 3rd order polynomial in cases where h has low slope, or where
% there is a local min/max in h
%
% TL: required a bit of thought to differentiate the polynomial regression
% step, see comments below.
tl_a_orig=tl_a;
warning off
for j=i_interp(:)'
  ind=unique(min(nx,max(1,j+[-2 -1 +1 +2])));
  ind=setdiff(ind,j);

  % pa=polyfit(ind,a_orig(ind),length(ind)-1);  % 3rd order polynomial
  % a(j)=polyval(pa,j);  % interpolate

  % Note: Here is the same thing as the above 2 lines of polyfit/polyval, but
  % written out in polynomial regression form.  This is verified to give the
  % same result, but writing this out this way makes it easier to add in the
  % TL code.
  for n=0:norder
    X(:,norder-n+1)=ind.^n;
  end
  Xi=pinv(X);  % note Xi is constant, not TL-dependent
  % pa = Xi*a_orig(ind)';
  tl_pa = Xi*tl_a_orig(ind)';
  % a(j)=0;
  tl_a(j)=0;
  for n=0:norder
    % a(j) = a(j) + pa(n+1)*j^(norder-n);
    tl_a(j) = tl_a(j) + tl_pa(n+1)*j^(norder-n);
  end

end
warning on

% interpolate to get q at predictor time step dt/2, eqn 20
tl_xp=-tl_a*dt/2;
tl_qp=nan(nx,1);
for i=1:nx
  i0=max(find(x<xp(i)));
  if(isempty(i0)) % extrapolate, x<x(1)
    x0=x(1)-dx;
    tl_q0=2*tl_q(2)-tl_q(1);
    tl_qp(i) = tl_q0 + (tl_q(1)-tl_q0)/dx*(xp(i)-x0) ...
        + (q(1)-q0)/dx*tl_xp(i);
  elseif(i0==nx)  % extrapolate, x>x(nx)
    x1=x(i0)+dx;
    tl_q1=2*tl_q(i0)-tl_q(i0-1);
    tl_qp(i) = tl_q(i0) + (tl_q1-tl_q(i0))/dx*(xp(i)-x(i0)) ...
        + (q1-q(i0))/dx*tl_xp(i);
  else  % interpolate
    tl_qp(i) = tl_q(i0) + (tl_q(i0+1)-tl_q(i0))/dx*(xp(i)-x(i0)) ...
        + (q(i0+1)-q(i0))/dx*tl_xp(i);
  end
end

% define candidate functions dh1,dh2,dh3, then apply limiters (MinMod
% function) to get dh, (h' in Marieu et al.), see non-numbered eqns
% following eqn 17
tl_dh1(2:nx) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
tl_dh1(1) = beta*(tl_h(2)-tl_h(1));
tl_dh2(2:nx-1) = .5*(tl_h(3:nx)-tl_h(1:nx-2));
tl_dh2(1) = tl_h(2)-tl_h(1);
tl_dh2(nx) = tl_h(nx)-tl_h(nx-1);
tl_dh3(1:nx-1) = beta*(tl_h(2:nx)-tl_h(1:nx-1));
tl_dh3(nx) = beta*(tl_h(nx)-tl_h(nx-1));
tl_dh=nan(nx,1);
for i=2:nx-1
  if((dh1(i)>0 & dh2(i)>0 & dh3(i)>0) || (dh1(i)<0 & dh2(i)<0 & dh3(i)<0))
    if(dh(i)==dh1(i))
      tl_dh(i)=tl_dh1(i);
    elseif(dh(i)==dh2(i))
      tl_dh(i)=tl_dh2(i);
    else
      tl_dh(i)=tl_dh3(i);
    end
  else
    tl_dh(i)=0;
  end
end
tl_dh(1)=tl_dh(2);
tl_dh(nx)=tl_dh(nx-1);

% h at time step t+dt, eqn 17.  Output is on the staggered grid
tl_hp=nan(nx-1,1);
tl_hp(1:nx-1) = .5*(tl_h(1:nx-1)+tl_h(2:nx)) ...
    + .125*(tl_dh(1:nx-1)-tl_dh(2:nx)) ...
    - dt/dx*(tl_qp(2:nx)-tl_qp(1:nx-1));

% output is supposed to be on the staggered grid, so interpolate q accordingly
tl_qpout=.5*(tl_qp(1:nx-1)+tl_qp(2:nx));

% revert the result to h=depth instead of h=elevation (see start of code
% where this was flipped)
tl_hpout=-tl_hp;
