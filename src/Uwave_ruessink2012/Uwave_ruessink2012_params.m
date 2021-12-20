function [Aw,Sw,Uw,workspc]=Uwave_ruessink2012_params(Hmo,k,omega,h)
%
% [Aw,Sw,Uw,workspc]=Uwave_ruessink2012_params(Hmo,k,omega,h)
%
% Helper function for Uwave_ruessink2012().  Calculates asymmetry (Aw),
% skewness (Sw), and amplitude (Uw).  This function can also be useful if
% you just want the parameters and don't want u(t).
%
% Hmo   : significant wave height, m
% k     : wavenumber, rad/m
% omega : wave frequency, rad/s
% h     : water depth, m
% 

% Ursell number, eqn (6)
aw=Hmo/2;
Ur=3/4*aw.*k./(k.*h).^3;

% wave velocity magnitude, linear theory, stated in text
Hrms=Hmo/1.4;
Uw = omega/2.*Hrms./sinh(k.*h);

% non-linearity param B, eqn (9).  Using parameter values quoted in text
p1=0;
p2= .857;  % +/- .016
p3=-.471;  % +/- .025
p4= .297;  % +/- .021
ee=(p3-log10(Ur))./p4;  % dec 17, 2021: This should use log10(), confirmed vs fig 2
dens=1+exp(ee);
B0 = p1 + (p2-p1)./dens;

% non-linearity phase param psi, eqn (10)
p5=.815;  % +/- .055
p6=.672;  % +/- .073
psi0 = -pi/2 + pi/2*tanh(p5./Ur.^p6);

% The above baseline B0 and psi0 give a baseline value for asymmetry and
% skewness (Aw and Sw)
Aw = B0.*sin(psi0);
Sw = B0.*cos(psi0);

% now save all relevant variables in a struct, so they can be reused in
% TL-AD functions
if(nargout>1)
  vname={};
  vname{end+1}='aw';
  vname{end+1}='Ur';
  vname{end+1}='Hrms';
  vname{end+1}='Uw ';
  vname{end+1}='p1';
  vname{end+1}='p2';
  vname{end+1}='p3';
  vname{end+1}='p4';
  vname{end+1}='p5';
  vname{end+1}='p6';
  vname{end+1}='ee';
  vname{end+1}='dens';
  vname{end+1}='B0';
  vname{end+1}='psi0';
  % vname{end+1}='B';
  % vname{end+1}='psi';
  % vname{end+1}='phi';
  % vname{end+1}='b';
  % vname{end+1}='r';
  % vname{end+1}='f';
  % vname{end+1}='f1';
  % vname{end+1}='f2';
  % vname{end+1}='u';
  vname{end+1}='Hmo';
  vname{end+1}='k';
  vname{end+1}='omega';
  vname{end+1}='h';
  vname{end+1}='Aw';
  vname{end+1}='Sw';
  workspc=struct;
  for i=1:length(vname)
    workspc=setfield(workspc,vname{i},eval(vname{i}));
  end
end
