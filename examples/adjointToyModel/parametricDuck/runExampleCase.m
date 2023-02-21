% Use RAH's bathymetry code to generate an example case.  Values are loosely
% based on the example of parametricBeaches_RAH/example1DCase.m, a Duck
% profile
%
addpath parametricBeaches_RAH
addpath(genpath('../../../src'))
clear

sedmodel='dubarbier';
% sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% define model grid
hmin=.25;
dx=5;
x=[100:dx:900]';
nx=length(x);

% generate bathymetry, with hard-coded inputs to make it "Duck-like"
xs=100;
betaShore=0.1;   % average shoreline beach slope (from literature)
xb=200;
xOff = 700; hOff = 7.5;     % chose a deep point from other information
betaOff = 0.0095;           % bathy slope at deep point
hSea = 4.5;                 % from Ruessink paper for Duck
h=make1DBeachEngine(x,xs,betaShore,xb,xOff,hOff,betaOff,hSea);
h(h<=hmin)=hmin;
h=flipud(h);

% model input parameters
d50=180e-6;
d90=400e-6;
H0=1.5;
theta0=0;
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

% run the NL model
[Hrms,vbar,theta,kabs,Q,hp,bkgd] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,...
                  windW,detady,d50,d90,params,sedmodel);
