function [hpout,qpout,bkgd]=dhdt_marieu2007(qin,hin,dx,dt,ghost)
%
% [hp,qp,bkgd]=dhdt_marieu2007(q,h,dx,dt,ghost)
%
% Modified nonoscillatory-centered scheme by Marieu et al. (2007, JGR)
% "Modeling of vortex ripple morphodynamics"
%
% Assumes x is monotonically increasing.
%
% x,q, and h should all be Nx1 row vectors, where N is the number of
% gridpoints.
%
% The code calculates updated bathymetry after time step dt, defined on a
% staggered grid with nx-1 gridpoints.  Thus, a ghost point needs to be
% added if you want nx output points.  The parameter 'ghost' should be
% given as +1, -1, or 0, with the following behaviour: 
%
%    ghost= 0: No ghost points are added, output has size nx-1
%    ghost=+1: Ghost point is added to end of grid, output point i corresponds to grid point i+1/2
%    ghost=-1: Ghost point is added to start of grid, output point i corresponds to grid point i-1/2
%
% Marieu et al. used a two-stage application with time step dt/2 to return
% back to the original (non-staggered) grid after time dt.  To acheive this,
% make two successive calls to this function with ghost=-1 for the first
% call and ghost=+1 for the second call, and dt/2 for each.  Pass the output
% h from the first call as the input to the second call.  Here is the
% recipe:
%
%   >> [hp1,qp1]=dhdt_marieu2007(q,h,dx,dt/2,-1);
%   >> hp=dhdt_marieu2007(qp1,hp1,dx,dt/2,+1);
%
% OPTIONAL: Output struct 'bkgd' contains background information required
% as input for TL-AD codes

% Add ghost points. Note here, Marieu's 'h' is elevation while mine is
% depth, so the sign of h is negated compared to eqn 17.  To make things
% simple, switch to their convention here then switch back to my convention
% at the end of the function
nxin=length(qin);
nx=nxin+1;
q=zeros(nx,1);
h=zeros(nx,1);
if(ghost==+1)
  q(1:nxin)=qin;
  q(nxin+1)=qin(nxin);
  h(1:nxin)=-hin;
  h(nxin+1)=-hin(nxin);
else
  q(1)=qin(1);
  q(2:(nxin+1))=qin;
  h(1)=-hin(1);
  h(2:(nxin+1))=-hin;
end
x=([1:nx]-1)*dx;

% flag potentially-bad points: where h has low slope, or where there is a
% local min/max in h
dhmin=0.01;
i_interp=find(abs(h(3:nx)-h(1:nx-2))<dhmin)+1;  % low slope
if(abs(h(2   )-h(1 ))<dhmin) i_interp=union(i_interp, 1); end  % low slope, boundary
if(abs(h(nx)-h(nx-1))<dhmin) i_interp=union(i_interp,nx); end  % low slope, boundary
i_interp = union(i_interp,...
                 find(sign(h(3:nx)-h(2:nx-1))==sign(h(1:nx-2)-h(2:nx-1)))+1);  % min/max
i_nointerp=setdiff(1:nx,i_interp);

% local bed celerity a, eqn (21), using centered difference for the valid
% points not flagged above
for i=i_nointerp(:)'
  if(i==1)
    a(i)=(q(i+1)-q(i))/(h(i+1)-h(i));
  elseif(i==nx)
    a(i)=(q(i)-q(i-1))/(h(i)-h(i-1));
  else
    a(i)=2*(q(i+1)-q(i-1))./(h(i+1)-h(i-1));
  end
end

% v1: Replace the potentially-bad points with 3rd order polynomial.  This
% version uses matlab's polyfit/polyval.
warning off
for j=i_interp(:)'
  ind=setdiff(1:nx,[j i_interp(:)']);
  thisdx=x(j)-x(ind);
  [~,isort]=sort(abs(thisdx),'ascend');
  ind=ind(isort(1:4));  % four nearest points
  pa=polyfit(ind,a(ind),length(ind)-1);  % 3rd order polynomial
  a(j)=polyval(pa,j);  % interpolate
end
warning on

% v2: Replace the potentially-bad points with 3rd order polynomial.  This
% version uses polynomial regression written out in matrix form.  This is
% verified to give the same result as v1 above (polyval/polyfit), but makes
% it easier to write the TL code.  Commented out of NL code after verifying
% that it matches the polyval/polyfit version.  Verification used duck94
% offshore bar migration case (Oct 4, 1994).
a2=a;
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
  pa = Xi*a(ind)';
  a2(j)=0;
  for n=0:norder
    a2(j) = a2(j) + pa(n+1)*j^(norder-n);
  end
end
if(max(abs(a2-a))>.01*max(abs(a)))
  error('mismatched polynomial regression solutions, should not happen')
end

% interpolate to get q at predictor time step dt/2, eqn 20
xp=x-a*dt/2;
qp=zeros(nx,1);
for i=1:nx
  i0=max(find(x<xp(i)));
  if(isempty(i0)) % extrapolate, x<x(1)
    qp(i)=q(1)+(q(2)-q(1))/dx*(xp(i)-x(1));
  elseif(i0==nx)  % extrapolate, x>x(nx)
    qp(i)=q(nx)+(q(nx)-q(nx-1))/dx*(xp(i)-x(nx));
  else  % interpolate
    qp(i)=q(i0)+(q(i0+1)-q(i0))/dx*(xp(i)-x(i0));
  end
end

% define candidate functions dh1,dh2,dh3, then apply limiters (MinMod
% function) to get dh, (h' in Marieu et al.), see non-numbered eqns
% following eqn 17
beta=4;
dh1(2:nx) = beta*(h(2:nx)-h(1:nx-1));
dh1(1) = beta*(h(2)-h(1));
dh2(2:nx-1) = .5*(h(3:nx)-h(1:nx-2));
dh2(1) = h(2)-h(1);
dh2(nx) = h(nx)-h(nx-1);
dh3(1:nx-1) = beta*(h(2:nx)-h(1:nx-1));
dh3(nx) = beta*(h(nx)-h(nx-1));
dh=zeros(nx,1);
for i=2:nx-1
  if(dh1(i)>0 & dh2(i)>0 & dh3(i)>0)
    dh(i)=min([dh1(i) dh2(i) dh3(i)]);
  elseif(dh1(i)<0 & dh2(i)<0 & dh3(i)<0)
    dh(i)=max([dh1(i) dh2(i) dh3(i)]);
  else
    dh(i)=0;
  end
end
dh(1)=dh(2);
dh(nx)=dh(nx-1);

% h at time step t+dt, eqn 17.  Output is on the staggered grid
hp=zeros(nx-1,1);
hp(1:nx-1) = .5*(h(1:nx-1)+h(2:nx)) ...
    + .125*(dh(1:nx-1)-dh(2:nx)) ...
    - dt/dx*(qp(2:nx)-qp(1:nx-1));

% output is supposed to be on the staggered grid, so interpolate q accordingly
qpout=.5*(qp(1:nx-1)+qp(2:nx));

% revert the result to h=depth instead of h=elevation (see start of code
% where this was flipped)
hpout=-hp;

% if requested, save info for use in TL-AD codes
if(nargout==3)
  bkgd=struct;
  bkgd.i_interp=i_interp;
  bkgd.i_nointerp=i_nointerp;
  bkgd.q       =q       ;
  bkgd.qin     =qin     ;
  bkgd.qp      =qp      ;
  bkgd.qpout   =qpout   ;
  bkgd.dh1     =dh1     ;
  bkgd.dh2     =dh2     ;
  bkgd.dh3     =dh3     ;
  bkgd.dh      =dh      ;
  bkgd.hp      =hp      ;
  bkgd.hpout   =hpout   ;
  bkgd.hin     =hin     ;
  bkgd.h       =h       ;
  bkgd.dx      =dx      ;
  bkgd.dt      =dt      ;
  bkgd.ghost   =ghost   ;
  bkgd.dhmin   =dhmin   ;
  bkgd.x       =x       ;
  bkgd.xp      =xp      ;
  bkgd.beta    =beta    ;
  bkgd.nx      =nx      ;
  bkgd.nxin    =nxin    ;
end
