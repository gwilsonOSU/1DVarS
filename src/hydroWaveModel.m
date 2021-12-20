function [Hrms,vbar,theta,kabs,Ew,Er,Dr,Aw,Sw,Uw,workspc] = ...
    hydroWaveModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,dAw,dSw)
%
% RECOMMENDED-USAGE: struct output
%
% out = hydroWaveModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,dAw,dSw)
%
% ALTERNATIVE-USAGE: individual variables output
%
% [Hrms,vbar,theta,kabs,Aw,Sw,Uw,workspc] = ...
%    hydroWaveModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,dAw,dSw)
%
% PURPOSE: Main front-end code for waves and currents.  If you also want
% morphology, use hydroSedModel.m.  This code is the same, but strips out
% the sediment transport.
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
% dAw      : linear correction factor for wave model Aw
% dSw      : linear correction factor for wave model Sw
% tau_wind : wind stress, vector Nx2, N/m2
% detady   : alongshore pressure gradient, m/m units
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
% Aw   : wave velocity asymmetry
% Sw   : wave velocity skewness
% Uw   : wave velocity amplitude

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

% wave shape parameters
[Aw0,Sw0,Uw,uwave_bkgd]=Uwave_ruessink2012_params(Hrms,kabs,omega,h);
Aw=Aw0+dAw;
Sw=Sw0+dSw;

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
vname{end+1}='dAw';  % input
vname{end+1}='dSw';  % input
vname{end+1}='imask';
vname{end+1}='Dr';
vname{end+1}='Er';
vname{end+1}='Ew';
vname{end+1}='Uw';
vname{end+1}='Aw';
vname{end+1}='Sw';
vname{end+1}='Hrms';
vname{end+1}='c';
vname{end+1}='h';
vname{end+1}='kabs';
vname{end+1}='nx';
vname{end+1}='theta';
vname{end+1}='hydro_bkgd';
vname{end+1}='uwave_bkgd';
workspc=struct;
for i=1:length(vname)
  if(exist(vname{i}))
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end

end  % main function for single time step
