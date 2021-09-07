function [Hrms,vbar,theta,kabs,Qx,hp,workspc] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,...
                  dgamma,...
                  d50,d90,params,sedmodel,dt)
%
% [Hrms,vbar,theta,kabs,Q,hp,wkspc]=hydroSedModel(x,h,H0,theta0,omega,ka_drag,...
%                                                 dgamma,...
%                                                 tau_wind,detady,d50,d90,sedparams,sedmodel,dt)
%
% Main front-end code for hydrodynamics and sediment transport
%
% INPUTS:
%
% x        : grid is +'ve onshore, and x(1) = offshore boundary
% h        : initial water depth, m
% H0       : rms wave height at offshore boundary, m
% theta0   : wave angle at offshore boundary, rads
% omega    : wave frequency, rad/m
% ka_drag  : hydraulic roughness factor, m
% dgamma   : linear correction factor for wave model gamma
% tau_wind : wind stress, vector Nx2, N/m2
% detady   : alongshore pressure gradient, m/m units
% d50      : median grain diameter, m
% d90      : 90th percentile grain diameter, m
% params.fv: breaking-induced eddy viscosity calibration parameter, see
%            Reniers et al. (2004) Table 4.  Scalar, recommended 0.1
% sedmodel : can be 'dubarbier', 'vanderA', or 'soulsbyVanRijn'.  See
%            corresponding codes qtrans_*.m and note their differing 'param'
%            struct, to be passed in here as part of 'params'
% dt       : time step (s) for integrating Exner equation to obtain update hp
%
% OUTPUTS:
%
% Hrms : rms wave height, m
% vbar : longshore current, m/s
% theta: wave angle, rads
% kabs : scalar wavenumber, rad/m
% Q    : cross-shore sediment flux, m2/s
% hp   : bathymetry updated after time step dt
% wkspc: all internal variables used in model; this is for passing as a
%        background state to TL-AD codes
%

% experimental features
doFilterQ=1;  % apply a filter to avoid sharp discontinuities in Q(x)
doDubarbierHack=0;  % shift velocities shoreward per Dubarbier's approxmiation
doMarieu=0;  % use Marieu's dh/dt instead of upwind differencing
if(doMarieu==1)
  warning('Marieu code TL-AD appears to have stability issues!')
end

physicalConstants;

horig=h;
imask=find(h<hmin);
h(imask)=hmin;  % min depth constraint

nx=length(x);

% % wind stress, following Reniers et al. (2004) eqn (6)
% cd=0.002;   % recommended in text
% w1=windW(:,1).^2+windW(:,2).^2;
% tau_wind(:,1)=cd*rhoa*sqrt(w1).*windW(:,1);
% tau_wind(:,2)=cd*rhoa*sqrt(w1).*windW(:,2);

% 1DH wave and longshore current balance
[Hrms,theta,vbar,kabs,Ew,Er,Dr,hydro_bkgd] = ...
    hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma);

% special case dt=0: in this case we are done since there is no morphology
% update to compute
if(dt==0)
  Qx=zeros(nx,1);
  hp=horig;
  h=horig;
else

% convert from (kabs,theta) to vector wavenumber, and calculate c for
% convenience
k=cat(2,kabs(:).*cos(theta(:)),kabs(:).*sin(theta(:)));
c=omega./kabs;

% depth averaged mean flow, Nx2 vector, +'ve onshore
ubar=-(Ew+2*Er)./(rho*c.*h);  % e.g., Dubarbier et al. (2015) eqn 8
ubar=cat(2,ubar(:),vbar(:));

% settling velocity: use Brown & Lawler.  TODO, for vanderA use 0.8*d50
% here, per explanation on page 29
if(strcmp(sedmodel,'vanderA'))
  ws=ws_brownLawler(.8*d50);
else
  ws=ws_brownLawler(d50);
end

% OPTIONAL: Dubarbier et al. suggest a modification to the mean velocity
% prior to calculation of undertow (udelta).  TODO, need TL-AD code
if(doDubarbierHack)
  warning('no TL-AD model exists for this optional code yet')
  lambda=1.57;
  xb=lambda*2*pi./kabs;
  for i=1:nx
    ind=find(x(i)-xb(i)<=x&x<=x(i));
    if(length(ind)<=1)
      ur(i,:)=ubar(i,:);
    else
      xx=xb(i)-x(i)+x(ind);
      ur(i,:)=trapz(x(ind),xx.*ubar(ind,:))/trapz(x(ind),xx);
    end
  end
  ubar=ur;
end

% Reniers et al. (2004) model for velocity at top of boundary layer
Dr(Dr==0)=min([0; Dr(Dr>0)]);
delta=zeros(nx,1);
for i=1:nx
  if(Dr(i)==0)
    udelta(i,:)=[0 0];
    delta(i)=0;
  else
    [udelta(i,:),delta(i),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),k(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           tau_wind(i,:),Dr(i),params.fv,d50(i));
  end
end

% rotate udelta into wave direction, as assumed by sed transport equations
udelta_w(:,1) = udelta(:,1).*cos(theta) - udelta(:,2).*sin(theta);
udelta_w(:,2) = udelta(:,1).*sin(theta) + udelta(:,2).*cos(theta);

% run the requested model for sediment flux (m2/s)
tanbeta=calcTanbeta(x,h)';
if(strcmp(sedmodel,'dubarbier'))  % Dubarbier et al. (2015)
  [Q,Qb,Qs,Qa,bkgd_qtrans] = ...
      qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta_w,ws,...
                       params.Cw,params.Cc,params.Cf,params.Ka);
elseif(strcmp(sedmodel,'soulsbyVanRijn'))  % Soulsby & van Rijn
  [Q,~] = ...
      qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,kabs,...
                            omega,theta,ubar,Dr,params);
elseif(strcmp(sedmodel,'vanderA'))  % van Der A et al. (2013)
  [Q,bkgd_qtrans] = ...
      qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta_w,delta,ws,params);
end
Q=real(Q);

% Dubarbier et al. suggest a threshold to exclude bed level change in the
% swash zone.  Use this, or the imposed minimum depth criteria, whichever
% is further offshore.
Trel=(2*pi/omega)*sqrt(g./h);
ishore=find(Trel>=40);
if(isempty(ishore))
  imax=nx;
else
  imax=min(ishore);
end
if(~isempty(imask))
  imax=min(min(imask),imax);
end
Q(imax:end)=0;

% OPTIONAL: use filtering to suppress instabilities.  This is necessary for
% vanderA model, to suppress discontinuities associated with switching
% on/off of phase lags, ripples, etc.
if(doFilterQ)
  fc = 1/50;   % band lower bound in 1/m
  fs=1/abs(x(2)-x(1));
  [b,a] = butter(5,fc/(fs/2),'low');
  Q = filtfilt(b,a,Q);
end

% TEST-CODE: weighting fn for damping out step function transport at shore.
% wgt=ones(nx,1);  % v1, no weighting
% wgt=(tanh((h-1.5)*2)+1)/2;  % v1, based on depth
% wgt=(1-tanh(([1:nx]'-min(imask))/2))/2;  % v2, based on shoreline
% wgt=(1-tanh(([1:nx]'-imax+2)/2))/2;  % v3, based on dubarbier condition
% Q=Q.*wgt;  % mitigate transport discontinuity at the shoreline

% weighting v4, cosine shape
% L=ceil(100/diff(x(1:2)));
L=20;
di=[1:nx]'-imax;
wgt=.5*(1+cos(pi/L*(di+L)));
wgt(di<-L)=1;
wgt(di>0)=0;

% rotate output from wave-following coords to cartesian
Qx=Q.*cos(theta);  % x-shore component

% bathymetry update: dhdt = -dzdt = dQdx.  This is the Exner equation,
% e.g. see Dubarbier et al. (2015) eqn. (16), and note Q is the volumetric
% transport rate (m2/s) including the bed porosity
if(doMarieu)  % use "stable" Marieu formulation for dh/dt
  dx=abs(x(2)-x(1));
  [hp1,qp1,bkgd_marieu_step1]=dhdt_marieu2007(Qx,horig,dx,dt/2,-1);
  [hp ,qp ,bkgd_marieu_step2]=dhdt_marieu2007(qp1,hp1,dx,dt/2,+1);
  dh=hp-horig;
  dQdx=nan(nx,1);
else
  dQdx=ddx_upwind(x,Qx,horig);
  dh=dQdx*dt;  % use ddx_upwind() result
  qp=nan(nx,1);
end
dh=dh.*wgt;   % apply damping near shore
dh(isnan(dh))=0;
hp = horig + dh;

end  % catch for special case dt==0

% save all relevant variables in a struct, so they can be reused in TL-AD
% functions
vname={};
vname{end+1}='x';       % input
vname{end+1}='h';       % input
vname{end+1}='H0';      % input
vname{end+1}='theta0';  % input
vname{end+1}='omega';   % input
vname{end+1}='ka_drag'; % input
vname{end+1}='tau_wind';% input
vname{end+1}='detady';  % input
vname{end+1}='dgamma';  % input
vname{end+1}='d50';     % input
vname{end+1}='d90';     % input
vname{end+1}='params';  % input
vname{end+1}='sedmodel';% input
vname{end+1}='dt';      % input
vname{end+1}='imask';
vname{end+1}='wgt';
vname{end+1}='horig';
vname{end+1}='Dr';
vname{end+1}='Er';
vname{end+1}='Ew';
vname{end+1}='Hrms';
vname{end+1}='Q';
vname{end+1}='Qx';
vname{end+1}='dQdx';
vname{end+1}='bkgd_qtrans';
vname{end+1}='c';
vname{end+1}='h';
vname{end+1}='hp';
vname{end+1}='dh';
vname{end+1}='k';
vname{end+1}='kabs';
vname{end+1}='nx';
vname{end+1}='tanbeta';
vname{end+1}='theta';
vname{end+1}='ubar';
vname{end+1}='udel_bkgd';
vname{end+1}='hydro_bkgd';
vname{end+1}='delta';
vname{end+1}='udelta';
vname{end+1}='udelta_w';
vname{end+1}='vbar';
vname{end+1}='ws';
vname{end+1}='doMarieu';
vname{end+1}='doFilterQ';
vname{end+1}='doDubarbierHack';
if(strcmp(sedmodel,'dubarbier'))
  vname{end+1}='Qb';
  vname{end+1}='Qs';
  vname{end+1}='Qa';
end
if(doMarieu)
  vname{end+1}='bkgd_marieu_step1';
  vname{end+1}='bkgd_marieu_step2';
end
workspc=struct;
for i=1:length(vname)
  if(exist(vname{i}))
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end

% optional, if user only wants the struct
if(nargout==1)
  Hrms=workspc;
end
