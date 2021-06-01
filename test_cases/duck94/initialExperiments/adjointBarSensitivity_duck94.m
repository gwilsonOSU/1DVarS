% Test cases based on canonical on/offshore migration cases from Duck94.
% The bathymetry are fitted to crab profiles by using Rob's parametric beach
% code.  The wave conditions are set equal to "nominal" values that occurred
% during the experiment, since I do not want to bother with the real
% time-dependent wave conditions.
%
addpath(genpath('~/work/unfunded_projects/sedimentTransport1D_TLAD/hydroSedModel/'))
addpath tools
addpath ../parametricDuck/parametricBeaches_RAH
clear

%--------------------------------------
% user inputs
%--------------------------------------

% duck94 case: set this to a,b,c, or d...  These follow Gallagher et
% al. (1998) four standard test cases.  All but case (b) are offshore bar
% migration events.
duck94Case='b';

% sediment model: uncomment one of the following...
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% model grid
hmin=.5;
dx=2.5;
x=[100:dx:900]';
nx=length(x);

%--------------------------------------
% Case-dependent settings.  Depends on input 'duck94Case'.  Do not edit.
%--------------------------------------

datadir='data/duck94_fulldataset';

% define:
%   dnum(1) = EST datenum for start of this case
%   dnum(2) = EST datenum for end of this case
switch duck94Case
 case 'a'
  dnum(1)=datenum(1994,9,1,19,0,0);
  dnum(2)=datenum(1994,9,5,19,0,0);
  beachParam.xs=145;
  beachParam.xb=210;
  beachParam.hSea=3;
 case 'b'
  dnum(1)=datenum(1994,9,22,19,0,0);
  dnum(2)=datenum(1994,9,27,19,0,0);
  beachParam.xs=145;
  beachParam.xb=270;
  beachParam.hSea=4
 case 'c'
  dnum(1)=datenum(1994,10,2,22,0,0);  % Gallagher et al. used 1900-2200, but
                                      % no shoreface sensors available then?
  dnum(2)=datenum(1994,10,4,16,0,0);
  beachParam.xs=145;
  beachParam.xb=235;
  beachParam.hSea=4;
 case 'd'
  dnum(1)=datenum(1994,10,10,22,0,0);
  dnum(2)=datenum(1994,10,15,22,0,0);
  beachParam.xs=145;
  beachParam.xb=265;
  beachParam.hSea=3.5
end

% some parametric beach profile params, don't need to vary from case to case
beachParam.betaShore=0.12;   % average shoreline beach slope (from literature)
beachParam.xOff = 730;
beachParam.hOff = 7.5;     % chose a deep point from other information
beachParam.betaOff = 0.0095;           % bathy slope at deep point

% read the .dof file for the case, to get sensor positions
i=1;  % set to 1 for initial profile, 2 for final profile
doffn=[datadir '/Dropbox/Duck94_doffs/' datestr(dnum(i),'mmddhhMM') '.dof'];
meafn=[datadir '/Dropbox/Duck94_mvar/'  datestr(dnum(i),'mmddhhMM') '.mea'];
[dof.s,dof.p,dof.u,dof.v,dof.t]=readDofs(doffn);
[mea.s,mea.p,mea.u,mea.v,mea.t,mea.meta]=readMea(meafn);

% only consider sonar .dof points for which there are existing .mea sensors
clear ind
for i=1:length(mea.s)
  for j=1:length(dof.s)
    if(strcmp(mea.s(i).name,dof.s(j).name))
      ind(i)=j;
      break
    end
  end
end
dof.s=dof.s(ind);

% take the 3-h average to get bathymetry points for this time
spuv.x=[mea.s.x];
spuv.zbed=mean(-[mea.s.data]+repmat([dof.s.zobs],[10 1]));

% average any non-unique points in sonar data
xnew=unique(spuv.x);
znew=nan*xnew;
for i=1:length(xnew)
  znew(i)=nanmean(spuv.zbed(spuv.x==xnew(i)));
end
spuv.x=xnew;
spuv.zbed=znew;
clear xnew znew

% TODO: define incident waves
warning('using hard-coded incident waves')
d50=180e-6;
d90=400e-6;
H0=.5;
omega=2*pi/8;
theta0=deg2rad(0);
ka_drag=0.015;
dt=24*60*60*7;

%--------------------------------------
% fit bathymetry to Holman et al. parametric shape
%--------------------------------------

beachParam0=beachParam;

% first cut, use hard-coded inputs that look close to what I want
h=make1DBeachEngine(x,...
                    beachParam.xs,...
                    beachParam.betaShore,...
                    beachParam.xb,...
                    beachParam.xOff,...
                    beachParam.hOff,...
                    beachParam.betaOff,...
                    beachParam.hSea);

% % refine parameters to fit data.  TODO: this seems to be unstable
% h=@(p)make1DBeachEngine(x,...
%                         p(1),... beachParam.xs,...
%                         p(2),... beachParam.betaShore,...
%                         p(3),... beachParam.xb,...
%                         beachParam.xOff,...
%                         beachParam.hOff,...
%                         beachParam.betaOff,...
%                         beachParam.hSea);
% hi=@(p)interp1(x,h(p),spuv.x);
% costFN=@(p)nansum( ( hi(p) - spuv.zbed ).^2 );
% pguess=[beachParam.xs;
%         beachParam.betaShore;
%         beachParam.xb];
% %         beachParam.xOff;
% %         beachParam.hOff;
% %         beachParam.betaOff;
% %         beachParam.hSea];
% p=fminsearch(costFN,pguess);
% beachParam.xs       =p(1);
% beachParam.betaShore=p(2);
% beachParam.xb       =p(3);
% % beachParam.xOff     =p(4);
% % beachParam.hOff     =p(5);
% % beachParam.betaOff  =p(6);
% % beachParam.hSea     =p(7);
% 
% % finalize the fitted parameterized bathymetry
% h=make1DBeachEngine(x,...
%                     beachParam.xs,...
%                     beachParam.betaShore,...
%                     beachParam.xb,...
%                     beachParam.xOff,...
%                     beachParam.hOff,...
%                     beachParam.betaOff,...
%                     beachParam.hSea);

% enforce minimum depth and handle NaN at shoreline
ishore=min(find(~isnan(h)));
h=h(ishore:end);
x=x(ishore:end);
h(h<hmin)=hmin;
nx=length(x);

% locate the bar extents x0 and x1, as the inflection point on either side
% of the profile
dhdx=diff(h);
dhdx(end+1)=dhdx(end);
dhdx2=diff(dhdx);
dhdx2(end+1)=dhdx2(end);
[~,ib]=min(abs(x-beachParam.xb));
i0=max(find([1:nx]'<ib & dhdx2<0));
i1=min(find([1:nx]'>ib & dhdx2<0));
x0=x(i0);
x1=x(i1);

% test plot: compare parameterized vs. observed bathy
clf, hold on
plot(x,-h)
plot(spuv.x,spuv.zbed,'r*-')
ax=axis;
plot([1 1]*x0,ax(3:4),'g--')
plot([1 1]*x1,ax(3:4),'g--')
box on

% finally, before proceeding note that the model expects x to be
% positive-onshore, so just flip the grid
h=flipud(h);
i0=nx-i0;
i1=nx-i1;

%--------------------------------------
% run the adjoint model
%--------------------------------------

% NL-model input parameters
tau_wind=zeros(nx,2);
detady=zeros(nx,1);
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
                  tau_wind,detady,d50,d90,params,sedmodel,dt);

% Adjoint sensitivity for bar migration.  "Cost function" = average sediment
% flux over the bar
J = .5*( Q(i0) + Q(i1) );
% tl_J = .5*( tl_Q(i0) + tl_Q(i1) );
ad_Q=zeros(nx,1);
ad_Q(i0)=.5;
ad_Q(i1)=.5;
ad_Hrms =zeros(nx,1);
ad_vbar =zeros(nx,1);
ad_theta=zeros(nx,1);
ad_kabs =zeros(nx,1);
ad_hp   =zeros(nx,1);
[ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,...
          ad_tau_wind,ad_detady,ad_d50,ad_d90,ad_params] = ...
    ad_hydroSedModel(ad_Hrms,ad_vbar,ad_theta,ad_kabs,ad_Q,ad_hp,bkgd);

% multiply adjoint variables to obtain sensitivities, based on 50-percent
% range around input values.  Exception is bathymetry 'feedback'
% sensitivity, for that I assume the variation is due to bathymetry change
% as predicted over the model time period dt (e.g., 1-hour)
disp(['sedmodel: ' sedmodel])
disp(['J: ' num2str(J                           ,'%e')])
disp(['feedback: ' num2str(sum(ad_h.*(hp-h))    ,'%e')])
disp(['H0     : ' num2str(ad_H0      *.5*H0     ,'%e')])
disp(['theta0 : ' num2str(ad_theta0  *.5*theta0 ,'%e')])
disp(['omega  : ' num2str(ad_omega   *.5*omega  ,'%e')])
disp(['ka_drag: ' num2str(ad_ka_drag *.5*ka_drag,'%e')])
disp(['d50    : ' num2str(ad_d50     *.5*d50    ,'%e')])
disp(['d90    : ' num2str(ad_d90     *.5*d90    ,'%e')])
for fld=fields(ad_params)'
  fld=cell2mat(fld);
  disp([fld ': ' num2str(getfield(ad_params,fld)*.5*getfield(params,fld),'%e')])
end

