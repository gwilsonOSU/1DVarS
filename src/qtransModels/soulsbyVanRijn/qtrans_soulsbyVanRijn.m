function [qtot,workspc] = qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,kabs,...
                                                omega,theta,ubar,Dr,param)
%
% [q,workspc] = qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,...
%                                              kabs,omega,theta,ubar,Dr,param)
%
% Calculates transport using the *basic* processes and parameterizations
% considered by xbeach, though with most of the fancy limiters and
% enhancements left out.  It is the Soulsby & van Rijn equilibrium sediment
% concentration, advected by a velocity representative of the mean flow and
% contributions from wave asymmetry & skewness.  Includes the simple
% "engineering" version of the bed slope effect also.
%
% Comments in code refer to equation numbers from the xbeach manual (March
% 2010 version), unless otherwise stated.
%
% INPUTS:
%
% All should be vector Nx1 unless otherwise stated
%
% x       : grid, m
% d50     : median grain size, m, scalar
% d90     : 90th percentile grain size, m, scalar
% h       : water depth, m
% tanbeta : bed slope
% Hrms    : rms wave height, m
% kabs    : scalar wavenumber, rad/m
% omega   : wave frequency, rad/m, scalar
% theta   : mean wave angle, radians
% ubar    : depth averaged below trough flow, Nx2 vector, m/s
% Dr      : roller dissipation, W/m2, used for diffusion
% param.alphab : calib. factor for slope-driven transport, Soulsby uses 1.6
% param.facua  : calib. factor for wave asymmetry, default 1.0
%

physicalConstants;

% drag coef for currents, from soulsby.  Note he says "z0 should be set to
% 6mm"
z0=.006;
Cd = (.4./(log(h/z0)-1)).^2;

% % old: advection velocity = mean flow, but add limiter for sheet flow, eqnq
% % 2.64.  In the end, I decided to omit this limiter to simplify things
% Delta=s-1;  % I think... they don't define it
% uE2_max =theta_sf*g*d50*Delta/cf;
% uE2_actual = sqrt( ubar(1)^2 + ubar(2)^2 );
% if(uE2_actual > uE2_max)
%   uE2 = uE2_max;
% else
%   uE2 = uE2_actual;
% end

% mean flow mag
uE2 = ubar(:,1).^2 + ubar(:,2).^2;

% critical velocity, Soulsby eqn 133d
ucr=0.19*d50^.1*log10(4*h./d90);

% rms wave orbital velocity.  Omits xbeach's finite amplitude "correction",
% and turbulence enhancement
urms = Hrms/2*omega./sinh(kabs.*h);

% wave asymmetry & skewness factors
Hmo = 1.4*Hrms;
aw=Hmo/2;
Ur=3/4*aw.*kabs./(kabs.*h).^3;
psi = -pi/2+pi/2*tanh(.64./Ur.^(.60));
earg = (-.61-log(Ur))/-.35;
dens = 1+exp(earg);
Sk = .79./dens.*cos(psi);
As = .79./dens.*sin(psi);

% advection velocity, plus extra term to represent wave asymmetry
VW = param.facua*urms.*(Sk-As);  % eqn 2.67
uAV = VW.*cos(theta)+ubar(:,1);  % eqn 2.66

% sediment flux in m2/s, eqns from Soulsby's text
Dstar=(g*(s-1)/nu^2)^(1/3)*d50;  % SCALAR
Asb = .005*h.*(d50./h).^1.2/((s-1)*g*d50)^1.2;  % eqn 136b
Ass = .012*d50*Dstar^(-.6)/((s-1)*g*d50)^1.2;  % eqn 136c, SCALAR
Ufact = sqrt(uE2 + .018*urms.^2./Cd) - ucr;
slopeFact = 1-param.alphab*tanbeta;
q=zeros(length(x),1);
C=zeros(length(x),1);
ind=find(Ufact>0);
q(ind) = (Ass+Asb(ind)).*uAV(ind).*Ufact(ind).^(2.4).*slopeFact(ind);
C(ind) = (Ass+Asb(ind)).*Ufact(ind).^(2.4).*slopeFact(ind)./h(ind);

% sediment diffusion coefficient Dh, this is copied from xbeach code since it's
% not specified in the manual nor the Roelvink paper
facDc=1;   % xbeach "turn on sediment diffusion" flag
nuhfac=1;  % xbeach "turn on roller-induced viscosity" flag
nuh=.1;    % bkgd viscosity, default 0.1 m2/s in xbeach
Dh = facDc*(nuh+nuhfac*h.*(Dr/rho).^(1/3));  % copied from xbeach code

% flux contribution from sediment diffusion
if(length(x)>2)
  dCdx=ddx_upwind(x,C);
  qtot = q + Dh.*dCdx;
else
  qtot=q;
end

% now save all relevant variables in a struct, so they can be reused in
% TL-AD functions
if(nargout>1)
  vname={};
  vname{end+1}='As';
  vname{end+1}='Asb';
  vname{end+1}='Ass';
  vname{end+1}='Cd';
  vname{end+1}='C';
  vname{end+1}='Dh';
  vname{end+1}='Dr';
  vname{end+1}='Dstar';
  vname{end+1}='Hmo';
  vname{end+1}='Hrms';
  vname{end+1}='Sk';
  vname{end+1}='Ufact';
  vname{end+1}='Ur';
  vname{end+1}='VW';
  vname{end+1}='aw';
  vname{end+1}='d50';
  vname{end+1}='d90';
  vname{end+1}='dens';
  vname{end+1}='earg';
  vname{end+1}='facDc';
  vname{end+1}='h';
  vname{end+1}='kabs';
  vname{end+1}='nuh';
  vname{end+1}='nuhfac';
  vname{end+1}='omega';
  vname{end+1}='param';
  vname{end+1}='psi';
  vname{end+1}='slopeFact';
  vname{end+1}='tanbeta';
  vname{end+1}='theta';
  vname{end+1}='uAV';
  vname{end+1}='uE2';
  vname{end+1}='ucr';
  vname{end+1}='ubar';
  vname{end+1}='urms';
  vname{end+1}='x';
  vname{end+1}='z0';
  vname{end+1}='q';
  workspc=struct;
  for i=1:length(vname)
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end
