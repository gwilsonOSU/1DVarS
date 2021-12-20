function [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Ew,tl_Er,tl_Dr,tl_Aw,tl_Sw,tl_Uw] = ...
    tl_hydroWaveModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,...
                     tl_detady,tl_dgamma,tl_dAw,tl_dSw,bkgd)%,outvar)


physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

nx=length(x);

%--------------------------------------
% begin TL
%--------------------------------------

tl_h(imask)=0;  % min depth constraint

% 1DH wave and longshore current balance
[tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
    tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,tl_detady,tl_dgamma,hydro_bkgd);

% wave shape parameters
[tl_Aw0,tl_Sw0,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hrms,tl_kabs,tl_omega,tl_h,uwave_bkgd);
tl_Aw=tl_Aw0+tl_dAw;
tl_Sw=tl_Sw0+tl_dSw;

% % TEST-CODE: override output variable
% if((length(outvar)>=6 & strcmp(outvar(1:6),'udelta')) | strcmp(outvar,'ubar'))
%   eval(['tl_Hrms = tl_' outvar '(:,2);']);
%   eval(['tl_' outvar '(:,1)=0;']);
% else
%   eval(['tl_Hrms = tl_' outvar ';']);
% end
% tl_vbar =zeros(nx,1);
% tl_theta=zeros(nx,1);
% tl_kabs =zeros(nx,1);
% tl_Qx   =zeros(nx,1);
% tl_hp   =zeros(nx,1);
