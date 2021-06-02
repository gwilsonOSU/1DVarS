% for a given beach profile, wave conditions, examine the sensitivity of
% on/offshore transport of a bar.  Requires knowing where the bar is.
%
addpath(genpath('../../src'))
clear

% background model
sedmodel='dubarbier';
% sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% define model grid
hmin=.25;
dx=5;
x=[100:dx:900]';
nx=length(x);

% generate bathymetry, with hard-coded inputs to make it "Duck-like"
addpath parametricDuck/parametricBeaches_RAH/
xs=100;
betaShore=0.1;   % average shoreline beach slope (from literature)
xb=200;
xOff = 700; hOff = 7.5;     % chose a deep point from other information
betaOff = 0.0095;           % bathy slope at deep point
hSea = 4.5;                 % from Ruessink paper for Duck
h=make1DBeachEngine(x,xs,betaShore,xb,xOff,hOff,betaOff,hSea);
h(h<=hmin)=hmin;
h=flipud(h);
plot(x,-h)
disp('click on bar extents')
xy=ginput(2);
x0=xy(1,1)
x1=xy(2,1)
[~,i0]=min(abs(x-x0));
[~,i1]=min(abs(x-x1));

% NL-model input parameters
d50=180e-6;
d90=400e-6;
H0=1.5;
theta0=deg2rad(10);
omega=2*pi/10;
windW=zeros(nx,2);
detady=zeros(nx,1);
ka_drag=0.015;
params.dt=60*60;
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
if(strcmp(sedmodel,'dubarbier'))
  params.Cw=0.00483 ;
  params.Cc=0.02002 ;
  params.Cf=0.01173 ;
  params.Ka=0.631e-4;
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;
  params.m=11;
  params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
  params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% bkgd solution from NL-model
[Hrms,vbar,theta,kabs,Q,hp,bkgd] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,...
                  windW,detady,d50,d90,params,sedmodel);

% Adjoint sensitivity for bar migration.  "Cost function" = average sediment
% flux over the bar
J = .5*( Q(i0) + Q(i1) );
% tl_J = .5*( tl_Q(i0) + tl_Q(i1) );
ad_Q=zeros(nx,1);
ad_Hrms =zeros(nx,1);
ad_vbar =zeros(nx,1);
ad_theta=zeros(nx,1);
ad_kabs =zeros(nx,1);
ad_hp   =zeros(nx,1);
ad_Q(i0)=.5;
ad_Q(i1)=.5;
[ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,...
          ad_windW,ad_detady,ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Q,ad_hp,bkgd);

