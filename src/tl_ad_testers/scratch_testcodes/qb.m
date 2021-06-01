function [q,workspc] = ...
    qb(tanbeta,h,Hrms,kabs,omega,udelta,Cw,Cc,Cf)

% other tunable parameters.  Here, hard-coded equal to the values used by
% Dubarbier et al. (2015)
eps_b=0.135;   % bedload sed efficiency param
tan_phi=0.63;  % sediment internal friction angle
nt=1000;  % number of timesteps per wave phase, used for calculating moments

% physical constants
physicalConstants;

% derived params
nx=length(h);
Hmo=Hrms*1.4;
Uw = omega.*Hrms./(2*sinh(kabs.*h));

% step through each gridpoint
t=linspace(0,2*pi/omega,nt);

for i=1:nx  % gridpoints
  [utilde(:,i),bkgd_uwave(i)]=Uwave_ruessink2012(omega*t,Hmo(i),kabs,omega,h(i));
  Aw(i) = omega*Uw(i);
  utot(:,i)=sqrt(utilde(:,i).^2+udelta(i,1)^2+udelta(i,2)^2);
  qb1(i)=mean( abs(utilde(:,i)).^2.*utilde(:,i) );
  qb2(i)=mean( abs(utilde(:,i)).^2*udelta(i,1) );
  qb3(i)=mean( abs(utot(:,i)).^2 );
  q(i) = rho*eps_b/tan_phi*( Cw*qb1(i) ...
                             + Cc*qb2(i) ...
                             - Cf*tanbeta(i)/tan_phi*qb3(i) );
end

% now save all relevant variables in a struct, so they can be reused in
% TL-AD functions
if(nargout>1)
  vname={};
  vname{end+1}='h';
  vname{end+1}='Hrms';
  vname{end+1}='omega';
  vname{end+1}='udelta';
  vname{end+1}='Cw';
  vname{end+1}='Cc';
  vname{end+1}='Cf';
  vname{end+1}='eps_b';
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
  vname{end+1}='Aw';
  vname{end+1}='utot';
  vname{end+1}='qb1';
  vname{end+1}='qb2';
  vname{end+1}='qb3';
  vname{end+1}='q';
  workspc=struct;
  for i=1:length(vname)
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end
