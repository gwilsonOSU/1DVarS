function [Q,Qb,Qs,Qa,workspc] = ...
    qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta,ws,...
                     Cw,Cc,Cf,Ka)
%
% [Q,Qb,Qs,Qa] = ...
%    qtrans_dubarbier(tanbeta,h,Hrms,kabs,omega,udelta,ws,[Cw,Cc,Cf,Ka])
%
% Implements the Unibest-like model outlined in Dubarbier et al. (2015)
% "Process-based modeling of cross-shore sandbar behavior"
%
% OUTPUTS: Q=Qb+Qs+Qa is total volumetric transport (m2/s units), udelta is
% the velocity (m/s) at top of boundary layer (Reniers et al., 2004).
%
% INPUTS: All are Nx1 vectors where N is the number of gridpoints, unless
% otherwise stated below
%
% tanbeta: Nx1 cross-shore grid of beach slope tan(beta)
% h     : water depth, m
% Hrms  : rms wave height, m
% kabs  : wavenumber, rad/m, scalar
% omega : wave frequency, rad/m, scalar
% udelta: near-bed mean flow (e.g., see udelta_reniers2004.m), Nx2 vector
%
% OPTIONAL: These tunable params will be given defaults if not provided
%
% Cw    : wave friction coef, scalar, default 0.00483 for Duck94
% Cc    : current friction coef, scalar, default 0.02002 for Duck94
% Cf    : gravitational friction coef, scalar, default 0.01173 for Duck94
% Ka    : acceleration skewness transport coef, scalar, default 0.631e-4 for Duck94
%
% Double checked all equations with Dubarbier et al. paper Apr 18, 2019

% other tunable parameters.  Here, hard-coded equal to the values used by
% Dubarbier et al. (2015), and these are the same used by Hsu et al. (2006)
eps_b=0.135;   % bedload sed efficiency param
eps_s=0.015;   % suspended sed efficiency param
tan_phi=0.63;  % sediment internal friction angle
nt=1000;  % number of timesteps per wave phase, used for calculating moments

% default values for input params
if(~exist('Cw')) Cw=0.00483;  end;
if(~exist('Cc')) Cc=0.02002;  end;
if(~exist('Cf')) Cf=0.01173;  end;
if(~exist('Ka')) Ka=0.631e-4; end;

% physical constants
physicalConstants;

% derived params
nx=length(h);
Hmo=Hrms*1.4;
Uw = omega.*Hrms./(2*sinh(kabs.*h));

% step through each gridpoint
t=linspace(0,2*pi/omega,nt);
for i=1:nx

  % intra-wave velocity information, using Ruessink et al. (2012).  Note,
  % comparing to Hsu et al. (2006) it seems as though Dubarbier et al. are
  % using \tilde{U} and \tilde{u} (upper and lower case) interchangably for
  % "wave velocity".  There are a few other obvious sloppy notational errors
  % like using "omega_s" for fall velocity
  [utilde(:,i),bkgd_uwave(i)]=Uwave_ruessink2012(omega*t,Hmo(i),kabs(i),omega,h(i));
  uH(:,i)=imag(hilbert(utilde(:,i)));
  Au_nums(i) = mean( uH(:,i).^3 );
  Au_dens_1(i) = mean( uH(:,i).^2 );  % == mean(utilde(:,i).^2)
  Au_dens(i) = Au_dens_1(:,i).^(3/2);
  Au(i) = Au_nums(i)/Au_dens(i);
  Aw(i) = omega*Uw(i);
  utot(:,i)=sqrt(utilde(:,i).^2+udelta(i,1)^2+udelta(i,2)^2);

  % energetics model, eqns (11)-(14).  There are a few very sloppy typos in
  % the paper, use Hsu et al. (2006) and Bailard (1981) for more coherent
  % references
  qa(i)=-(rhop-rho)*g*Ka*Au(i)*Aw(i);
  qb1(i)=mean( abs(utilde(:,i)).^2.*utilde(:,i) );
  qb2(i)=mean( abs(utot(:,i)).^2*udelta(i,1) );
  qb3(i)=mean( abs(utot(:,i)).^3 );
  qb(i) = rho*eps_b/tan_phi*( Cw*qb1(i) ...
                              + Cc*qb2(i) ...
                              - Cf*tanbeta(i)/tan_phi*qb3(i) );
  qs1(i)=mean(abs(utilde(:,i)).^3.*utilde(:,i));
  qs2(i)=mean(abs(utot(:,i)).^3*udelta(i,1));
  qs3(i)=mean(abs(utot(:,i)).^5);
  qs(i) = rho*eps_s/ws(i)*( Cw*qs1(i) ...
                         + Cc*qs2(i) ...
                         - Cf*eps_s*tanbeta(i)/ws(i)*qs3(i) );
  q(i)=qa(i)+qb(i)+qs(i);

end

% normalize to get volumetric transport, m2/s units, and convert to row
% vectors
qnorm=g*(rhop-rho)*(1-psed);
Qa=qa(:)/qnorm;
Qb=qb(:)/qnorm;
Qs=qs(:)/qnorm;
Q =q (:)/qnorm;

% now save all relevant variables in a struct, so they can be reused in
% TL-AD functions
if(nargout>4)
  vname={};
  vname{end+1}='h';
  vname{end+1}='Hrms';
  vname{end+1}='omega';
  vname{end+1}='udelta';
  vname{end+1}='ws';
  vname{end+1}='Cw';
  vname{end+1}='Cc';
  vname{end+1}='Cf';
  vname{end+1}='Ka';
  vname{end+1}='eps_b';
  vname{end+1}='eps_s';
  vname{end+1}='tan_phi';
  vname{end+1}='nt';
  vname{end+1}='nx';
  vname{end+1}='Hmo';
  vname{end+1}='kabs';
  vname{end+1}='Uw ';
  vname{end+1}='tanbeta';
  vname{end+1}='t';
  vname{end+1}='utilde';
  vname{end+1}='bkgd_uwave';
  vname{end+1}='uH';
  vname{end+1}='Au_nums';
  vname{end+1}='Au_dens_1';
  vname{end+1}='Au_dens';
  vname{end+1}='Au';
  vname{end+1}='Aw';
  vname{end+1}='utot';
  vname{end+1}='qb1';
  vname{end+1}='qb2';
  vname{end+1}='qb3';
  vname{end+1}='qs1';
  vname{end+1}='qs2';
  vname{end+1}='qs3';
  vname{end+1}='qa';
  vname{end+1}='qs';
  vname{end+1}='qb';
  vname{end+1}='qnorm';
  vname{end+1}='Qa';
  vname{end+1}='Qs';
  vname{end+1}='Qb';
  workspc=struct;
  for i=1:length(vname)
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end
