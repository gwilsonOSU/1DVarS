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

% Note, Marieu's 'h' is elevation while mine is depth, so the sign of h is
% negated compared to eqn 17.  To make things simple, switch to their
% convention here then switch back to my convention at the end of the
% function
h=-hin;

% add ghost points
nxin=length(qin);
if(ghost==+1)
  q=[qin; qin(nxin)];
  h=[hin; hin(nxin)];
else
  q=[qin(1); qin];
  h=[hin(1); hin];
end

% grid
nx=length(q);
x=([1:nx]-1)*dx;

% local bed celerity a, eqn (21), using centered difference
a(2:nx-1)=2*(q(3:nx)-q(1:nx-2))./(h(3:nx)-h(1:nx-2));
a(1)=(q(2)-q(1))/(h(2)-h(1));
a(nx)=(q(nx)-q(nx-1))/(h(nx)-h(nx-1));

% replace with 3rd order polynomial in cases where h has low slope, or where
% there is a local min/max in h
a_orig=a;
dhmin=0.01;
i_interp=find(abs(h(3:nx)-h(1:nx-2))<dhmin)+1;  % low slope
if(abs(h(2   )-h(1 ))<dhmin) i_interp=union(i_interp,1); end  % low slope, boundary
if(abs(h(nx)-h(nx-1))<dhmin) i_interp=union(nx,1); end  % low slope, boundary
i_interp=union(i_interp,find(sign(h(3:nx)-h(2:nx-1))==sign(h(1:nx-2)-h(2:nx-1)))+1);  % min/max
warning off
for j=i_interp(:)'
  ind=unique(min(nx,max(1,j+[-2 -1 +1 +2])));
  ind=setdiff(ind,j);
  pa=polyfit(ind,a_orig(ind),length(ind)-1);  % 3rd order polynomial
  a(j)=polyval(pa,j);  % interpolate
end
warning on

% interpolate to get q at predictor time step dt/2, eqn 20
xp=x-a*dt/2;
qp=nan(nx,1);
for i=1:nx
  i0=max(find(x<xp(i)));
  if(isempty(i0)) % extrapolate, x<x(1)
    x0=x(1)-dx;
    q0=2*q(2)-q(1);
    qp(i)=q0+(q(1)-q0)/dx*(xp(i)-x0);
  elseif(i0==nx)  % extrapolate, x>x(nx)
    x1=x(i0)+dx;
    q1=2*q(i0)-q(i0-1);
    qp(i)=q(i0)+(q1-q(i0))/dx*(xp(i)-x(i0));
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
dh=nan(nx,1);
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
hp=nan(nx-1,1);
hp(1:nx-1) = .5*(h(1:nx-1)+h(2:nx)) ...
    + .125*(dh(1:nx-1)-dh(2:nx)) ...
    - dt/dx*(qp(2:nx)-qp(1:nx-1));

% output is supposed to be on the staggered grid, so interpolate q accordingly
qpout=.5*(qp(1:nx-1)+qp(2:nx));

% % OLD: interpolate to non-staggered grid
% hpp=nan(nx,1);
% hpp(2:nx-1) = .5*(hp(1:nx-2)+hp(2:nx-1));
% hpp(1) = 1.5*hp(1) - .5*hp(2);
% hpp(nx) = 1.5*hp(nx) - .5*hp(nx-1);

% revert the result to h=depth instead of h=elevation (see start of code
% where this was flipped)
hpout=-hp;

% if requested, save info for use in TL-AD codes
if(nargout==3)
  bkgd=struct;
  bkgd.a_orig  =a_orig  ;
  bkgd.i_interp=i_interp;
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
