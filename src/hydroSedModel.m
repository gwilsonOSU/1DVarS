function [Hrms,vbar,theta,kabs,Qx,hpout,workspc] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
                  d50,d90,params,sedmodel,gammaType,betaType,dt,nsubsteps)%,outvar)
%
% RECOMMENDED-USAGE: struct output
%
% out =     hydroSedModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
%                         d50,d90,params,sedmodel,gammaType,betaType,dt,nsubsteps)
%
%
% ALTERNATIVE-USAGE: individual variables output
%
% [Hrms,vbar,theta,kabs,Q,hp,out] = hydroSedModel(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
%                                                 d50,d90,sedparams,sedmodel,gammaType,betaType,dt,nsubsteps)
%
%
% PURPOSE: Main front-end code for hydrodynamics and sediment transport
%
% INPUTS:
%
% x        : grid is +'ve onshore, and x(1) = offshore boundary
% h        : initial water depth, m
% H0       : rms wave height at offshore boundary, m
% theta0   : wave angle at offshore boundary, rads
% omega    : wave frequency, rad/m
% ka_drag  : hydraulic roughness factor, m
% beta0    : assuming betaType=='const', sets the value of beta for roller model
% tau_wind : wind stress, vector Nx2, N/m2
% detady   : alongshore pressure gradient, m/m units
% dgamma   : linear correction factor for wave model gamma
% dAw      : linear correction factor for wave model Aw
% dSw      : linear correction factor for wave model Sw
% d50      : median grain diameter, m
% d90      : 90th percentile grain diameter, m
% params.fv: breaking-induced eddy viscosity calibration parameter, see
%            Reniers et al. (2004) Table 4.  Scalar, recommended 0.1
% params.lambda: Dubarbier et al.'s eqn 9 to shift velocities shoreward.
%                Set lambda=0 to turn this off.  Default lambda=1.57.
% params.* : other parameters as required by sediment transport model (see
%            below)
% sedmodel : can be 'dubarbier', 'vanderA', or 'soulsbyVanRijn'.  See
%            corresponding codes qtrans_*.m and note their differing 'param'
%            struct, to be passed in here as part of 'params'
% gammaType: can be 2001 or 2003 (numbers), for wave breaking gamma
%            formulation used by Ruessink et al. 2001 or 2003
% betaType : sets roller model beta(x) formulation.  Can be 'const', 'none',
%             or 'rafati21' (variable beta model of Rafati et al., 2001)
% dt       : time step (s) for integrating Exner equation to obtain update hp
% nsubsteps: number of sub-steps to divide integration over dt, for
%            numerical stability.  Each substep uses the same boundary
%            conditions and input params, but bathymetry and hydrodynamcs
%            are updated over time
%
% OUTPUTS:
%
% out: all internal variables used in model; this is for passing as a
%      background state to TL-AD codes.  The variables listed below are
%      included in this struct.  If you use nsubsteps>1, then 'out' will be
%      an array of structs, one for each time step.
% Hrms : rms wave height, m
% vbar : longshore current, m/s
% theta: wave angle, rads
% kabs : scalar wavenumber, rad/m
% Q    : cross-shore sediment flux, m2/s
% hp   : bathymetry for each sub-timestep (0,1,...,nsubsteps)
%

if(~exist('nsubsteps'))
  nsubsteps=1;
end

% sub-stepping loop
hp(:,1) = h;  % init t=0
for n=1:nsubsteps
  % disp(['hydroSedModel substep ' num2str(n) ' of ' num2str(nsubsteps)])
  [Hrms(:,n),vbar(:,n),theta(:,n),kabs(:,n),...
   Qx(:,n),hp(:,n+1),workspc(n)] = ...
      hydroSedModel_main(x,hp(:,n),H0,theta0,omega,ka_drag,beta0,tau_wind,detady,...
                         dgamma,dAw,dSw,d50,d90,params,sedmodel,gammaType,betaType,dt/nsubsteps);%,outvar);
end
for n=1:nsubsteps
  workspc(n).nsubsteps=nsubsteps;
end

% drop initial condition from hpout
for n=1:nsubsteps
  hpout(:,n)=hp(:,n+1);
end

% handle output format: user may request to get just the struct
% (nargout==1), or they may want individual variables (nargout>1).
if(nargout==1)
  Hrms=workspc;
end

end  % end wrapper function (for sub-stepping loop logic)

% begin main function, for single time step
function [Hrms,vbar,theta,kabs,Qx,hp,workspc] = ...
    hydroSedModel_main(x,h,H0,theta0,omega,ka_drag,beta0,tau_wind,detady,dgamma,dAw,dSw,...
                  d50,d90,params,sedmodel,gammaType,betaType,dt)%,outvar)

% experimental features
doFilterQ=0;  % numerical filtering with low-pass filter.  Don't use unless nuQ==0 (either choose one filtering scheme, or turn both off)
nuQ=.1;  % horizontal diffusion coef to provide numerical filtering for Q(x).  Set to zero to turn off filtering.
nuN=50; % increase this number to increase degree of horizontal diffusive filtering
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
    hydro_ruessink2001(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,beta0,gammaType,betaType);

% don't allow Dr==0, otherwise the undertow model will be discontinuous at
% the break point
Drmin=.001;
Drmask=find(Dr<Drmin);
Dr(Drmask)=.001;

% apply masking to hydro outputs before proceeding to sediment transport
Hrms (imask)=0;
vbar (imask)=0;
Ew   (imask)=0;
Er   (imask)=0;
Dr   (imask)=0;

% wave shape parameters.  Note Uwave_ruessink2012 specifies Hmo as input
Hmo=1.4*Hrms;
[Aw0,Sw0,Uw,uwave_bkgd]=Uwave_ruessink2012_params(Hmo,kabs,omega,h);
Aw=Aw0+dAw;
Sw=Sw0+dSw;

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
ubarx=-(Ew+2*Er)./(rho*c.*h);  % e.g., Dubarbier et al. (2015) eqn 8
ubar0(:,1)=ubarx;
ubar0(:,2)=vbar;

% settling velocity: Brown & Lawler.  For vanderA use 0.8*d50 per
% explanation on page 29
if(strcmp(sedmodel,'vanderA'))
  ws=ws_brownLawler(.8*d50);
else
  ws=ws_brownLawler(d50);
end

% OPTIONAL: Dubarbier et al. suggest a modification to the mean velocity
% prior to calculation of undertow (udelta)
if(params.lambda>0)
  ubar=dubarbierUmod(ubar0,kabs,params.lambda,x);
else
  ubar=ubar0;
end

% Reniers et al. (2004) model for velocity at top of boundary layer
delta_bl=ones(nx,1)*.05;  % init
udelta=ubar;   % init
for i=1:nx
  if(Dr(i)>0)
    [udelta(i,:),delta_bl(i),udel_bkgd(i)]= ...
        udelta_reniers2004(ubar(i,:),k(i,:),omega,...
                           h(i),Hrms(i),detady(i),...
                           tau_wind(i,:),Dr(i),params.fv,params.ks,d50(i));
    udel_bkgd=udel_bkgd(:);
  end
end

% Rotate udelta into wave direction, as assumed by sed transport equations.
% To understand which direction to rotate, consider the first component of
% udelta_w(:,1) is the "cross-shore" current in the wave-rotated
% coordinates, and it should equal the projection of udelta onto the
% wave-direction unit vector (cos(theta),sin(theta)), which establishes the
% eqn for udelta_w(:,1).  The eqn for udelta_w(:,2) then follows because
% it has to be a rotation matrix.
udelta_w(:,1) = +udelta(:,1).*cos(theta) + udelta(:,2).*sin(theta);
udelta_w(:,2) = -udelta(:,1).*sin(theta) + udelta(:,2).*cos(theta);

% run the requested model for sediment flux (m2/s)
tanbeta=calcTanbeta(x,h)';
if(strcmp(sedmodel,'dubarbier'))  % Dubarbier et al. (2015)
  [Q0,Qb,Qs,Qa,bkgd_qtrans] = ...
      qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta_w,ws,Aw,Sw,Uw,...
                       params.Cw,params.Cc,params.Cf,params.Ka);
elseif(strcmp(sedmodel,'soulsbyVanRijn'))  % Soulsby & van Rijn
  [Q0,~] = ...
      qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,kabs,...
                            omega,theta,ubar,Dr,params);
elseif(strcmp(sedmodel,'vanderA'))  % van Der A et al. (2013)
  [Q0,bkgd_qtrans] = ...
      qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta_w,delta_bl,ws,Aw,Sw,Uw,params);
end
Q0=real(Q0);

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
imax=imax-1;
Q1=Q0;
Q1(imax:end)=0;

% OPTIONAL: use filtering to suppress instabilities.  This is necessary for
% vanderA model, to suppress discontinuities associated with switching
% on/off of phase lags, ripples, etc.
if(doFilterQ)
  if(nuQ>0)
    warning('Recommended to NOT use both doFilterQ and nuQ>0, instead choose one or the other')
  end
  fc = 1/50;   % band lower bound in 1/m
  fs=1/abs(x(2)-x(1));
  [b,a] = butter(5,fc/(fs/2),'low');
  Q=zeros(size(Q1));
  Q(1:(imax-1)) = filtfilt(b,a,Q1(1:(imax-1)));
else
  Q = Q1;
end

% OPTIONAL: Apply horizontal diffusion to Q(x)
if(nuQ>0)
  dx=diff(x(1:2));
  A=zeros(nx);
  for i=2:nx-1
    A(i,i+[-1:1])=[1 -2 1]/dx^2*nuQ;
  end
  A(1,1:2)=0;%[-2 1]/dx^2*nuQ;
  A(nx,nx-1:nx)=[1 -2]/dx^2*nuQ;
  diffusionSmoother=(eye(nx)+A)^nuN;
  Q = diffusionSmoother*Q1;
else
  Q=Q1;
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

% rotate output from wave-following coords back to cartesian
Qx=Q.*cos(theta);  % x-shore component

% apply masking
Qx(imask)=0;

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
  % dQdx=ddx_upwind(x,Qx,horig);
  dQdx=ddx_centered(x,Qx);
  dh=dQdx*dt;
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
vname{end+1}='beta0';   % input
vname{end+1}='tau_wind';% input
vname{end+1}='detady';  % input
vname{end+1}='dgamma';  % input
vname{end+1}='dAw';  % input
vname{end+1}='dSw';  % input
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
vname{end+1}='Uw';
vname{end+1}='Aw';
vname{end+1}='Sw';
vname{end+1}='Hrms';
vname{end+1}='Q';
vname{end+1}='Q0';
vname{end+1}='Q1';
vname{end+1}='imax';
vname{end+1}='Qx';
vname{end+1}='dQdx';
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
vname{end+1}='ubar0';
vname{end+1}='udel_bkgd';
vname{end+1}='hydro_bkgd';
vname{end+1}='uwave_bkgd';
vname{end+1}='bkgd_qtrans';
vname{end+1}='delta_bl';
vname{end+1}='udelta';
vname{end+1}='udelta_w';
vname{end+1}='vbar';
vname{end+1}='ws';
vname{end+1}='doMarieu';
vname{end+1}='doFilterQ';
vname{end+1}='nuQ';
vname{end+1}='nuN';
vname{end+1}='Drmin';
vname{end+1}='Drmask';
if(nuQ>0)
  vname{end+1}='diffusionSmoother';
end
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

% TEST: tweak output for TL testing
% eval(['Qx=' outvar ';']);

end  % main function for single time step
