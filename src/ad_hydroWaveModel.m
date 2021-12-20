function [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tau_wind,...
          ad_detady,ad_dgamma,ad_dAw,ad_dSw] = ...
    ad_hydroWaveModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Ew,ad_Er,ad_Dr,ad_Aw,ad_Sw,ad_Uw,bkgd)%,invar)

physicalConstants;

% get NL variables from saved workspace struct
fld=fields(bkgd);
for i=1:length(fld)
  eval([fld{i} ' = bkgd.' fld{i} ';']);
end

nx=length(h);

%----------------------------------------
% init AD
%----------------------------------------

ad_h=zeros(nx,1);
ad_H0=0;
ad_theta0=0;
ad_omega=0;
ad_ka_drag=0;
ad_tau_wind=zeros(nx,2);
ad_detady=zeros(nx,1);
ad_dgamma=zeros(nx,1);
ad_dAw=zeros(nx,1);
ad_dSw=zeros(nx,1);
ad_Aw0=zeros(nx,1);
ad_Sw0=zeros(nx,1);

%----------------------------------------
% begin AD code
%----------------------------------------

% wave shape parameters
% tl_Sw=tl_Sw0+tl_dSw;
ad_Sw0=ad_Sw0+ad_Sw;
ad_dSw=ad_dSw+ad_Sw;
ad_Sw=zeros(nx,1);
% tl_Aw=tl_Aw0+tl_dAw;
ad_Aw0=ad_Aw0+ad_Aw;
ad_dAw=ad_dAw+ad_Aw;
ad_Aw=zeros(nx,1);
% [tl_Aw0,tl_Sw0,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hrms,tl_kabs,tl_omega,tl_h,uwave_bkgd);
[ad1_Hrms,ad1_kabs,ad1_omega,ad1_h]=ad_Uwave_ruessink2012_params(ad_Aw0,ad_Sw0,ad_Uw,uwave_bkgd);
ad_Hrms =ad_Hrms +ad1_Hrms ;
ad_kabs =ad_kabs +ad1_kabs ;
ad_omega=ad_omega+ad1_omega;
ad_h    =ad_h    +ad1_h    ;

% 1DH wave and longshore current balance
% [tl_Hrms,tl_theta,tl_vbar,tl_kabs,tl_Ew,tl_Er,tl_Dr] = ...
%     tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,tl_tau_wind,tl_detady,tl_dgamma,hydro_bkgd);
[ad1_h,ad1_H0,ad1_theta0,ad1_omega,ad1_ka_drag,ad1_tau_wind,ad1_detady,ad1_dgamma] = ...
    ad_hydro_ruessink2001(ad_Hrms,ad_theta,ad_vbar,ad_kabs,ad_Ew,ad_Er,ad_Dr,hydro_bkgd);
ad_h       =ad_h       +ad1_h       ;
ad_H0      =ad_H0      +ad1_H0      ;
ad_theta0  =ad_theta0  +ad1_theta0  ;
ad_omega   =ad_omega   +ad1_omega   ;
ad_ka_drag =ad_ka_drag +ad1_ka_drag ;
ad_tau_wind=ad_tau_wind+ad1_tau_wind;
ad_detady  =ad_detady  +ad1_detady  ;
ad_dgamma  =ad_dgamma  +ad1_dgamma  ;

%1 tl_h(imask)=0;  % min depth constraint
ad_h(imask)=0;
