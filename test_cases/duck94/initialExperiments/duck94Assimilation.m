% Test cases of sequential data assimilation during canonical on/offshore
% migration cases from Duck94.
%
addpath(genpath('~/work/unfunded_projects/sedimentTransport1D_TLAD/hydroSedModel/'))
addpath tools
clear

%--------------------------------------
% user inputs
%--------------------------------------

% duck94 case: set this to a,b,c, or d...  These follow Gallagher et
% al. (1998) four standard test cases.  All but case (b) are offshore bar
% migration events.

duck94Case='b';

% sediment model: uncomment one of the following...
sedmodel='dubarbier';
% sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

% model grid
grid.dx=2.5;
grid.x=[100:grid.dx:900]';
grid.nx=length(grid.x);

%--------------------------------------
% Case-dependent settings.  Depends on input 'duck94Case'.  Do not edit.
%--------------------------------------

datadir='data/duck94_fulldataset';

% define:
%   dnum(1) = EST datenum for start of this case
%   dnum(2) = EST datenum for end of this case
%   bathyfn = full-domain CRAB survey, used for initial bathymetry offshore of SPUV-array
switch duck94Case
 case 'a'
  dnum(1)=datenum(1994,9,1,19,0,0);
  dnum(2)=datenum(1994,9,5,19,0,0);
  % bathyfn='data/crab/FRF_19940907_0705_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0907_profile';
 case 'b'
  dnum(1)=datenum(1994,9,22,19,0,0);
  dnum(2)=datenum(1994,9,27,19,0,0);
  % bathyfn='data/crab/FRF_19940921_0710_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0921_profile';
 case 'c'
  dnum(1)=datenum(1994,10,2,22,0,0);  % Gallagher et al. used 1900-2200, but
                                      % no shoreface sensors available then?
  dnum(2)=datenum(1994,10,4,16,0,0);
  % warning('test code, changed start time'); dnum(1)=728570.3;  % TEST CODE, start at a time with non-trivial waves
  % bathyfn='data/crab/FRF_19941004_0716_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1003_profile';
 case 'd'
  dnum(1)=datenum(1994,10,10,22,0,0);
  dnum(2)=datenum(1994,10,15,22,0,0);
  % bathyfn='data/crab/FRF_19941010_0722_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1021_profile';
end

%--------------------------------------
% get SPUV measurements for this time period
%--------------------------------------

% % load all .mea files to get SPUV data vs. time
% meafn=fileList_ls([datadir '/Dropbox/Duck94_mvar/{09,10}*.mea']);
% for i=1:length(meafn)
%   disp(['reading mea file ' num2str(i) ' of ' num2str(length(meafn))])
%   clear this
%   [this.s,this.p,this.u,this.v,this.t,this.meta]=readMea(meafn{i});
%   this.fn=meafn{i};
%   this.dnum_est=this.meta.dnum_est;
%   mea(i)=this;
% end
% save data/allmea.mat mea
load data/allmea.mat

% % load all .dof files, needed to convert sonar data to bottom elevation
% doffn=fileList_ls([datadir '/Dropbox/Duck94_doffs/{09,10}*.dof']);
% for i=1:length(doffn)
%   disp(['reading dof file ' num2str(i) ' of ' num2str(length(doffn))])
%   clear this
%   [this.s,this.p,this.u,this.v,this.t]=readDofs(doffn{i});
%   this.fn=doffn{i};
%   dstr=strsh(doffn{i},'tr');
%   this.dnum_est=datenum(['1994-' dstr(1:2) '-' dstr(3:4) ' ' ...
%                       dstr(5:6) ':' dstr(7:8)]);
%   dof(i)=this;
% end
% save data/alldof.mat dof
load data/alldof.mat

% % load all .hsg files to get Hsig vs. time
% hsgfn=fileList_ls([datadir '/Dropbox/Duck94_mvar/{09,10}*.hsg']);
% for i=1:length(hsgfn)
%   disp(['reading hsg file ' num2str(i) ' of ' num2str(length(hsgfn))])
%   clear this
%   [this.Hs,this.meta]=readHsg(hsgfn{i});
%   this.fn=hsgfn{i};
%   this.dnum_est=this.meta.dnum_est;
%   hsg(i)=this;
% end
% save data/allhsg.mat hsg
load data/allhsg.mat

% match up the dof timestamps with the mea ones
[~,i,j]=intersect([dof.dnum_est],[mea.dnum_est]);
dof=dof(i);

% the hsg and mea files all have matching time stamps, so can just merge the data together
if(max(abs([mea.dnum_est]-[hsg.dnum_est]))>0 | max(abs([mea.dnum_est]-[dof.dnum_est]))>0)
  error('should not happen')
end
for i=1:length(hsg)
  mea(i).H=hsg(i).Hs;
  for j=1:length(mea(i).H)
    mea(i).H(j).data=mea(i).H(j).data/1.4;  % convert Hsig to Hrms
  end
end
clear hsg

% go through the dofs and transfer the relevant data into the mea struct.
% This is needed for calculating elevations from sonar altimeter data.  Note
% the timestamps are already aligned so that mea(i) corresponds to dof(i)
for i=1:length(mea)
  for vname='spuvt'
    thismea=getfield(mea(i),vname);
    thisdof=getfield(dof(i),vname);
    for j=1:length(thismea)
      for k=1:length(thisdof)
        if(strcmp(thismea(j).name,thisdof(k).name))
          break;
        end
      end
      thismea(j).zobs=thisdof(k).zobs;
      thismea(j).zbed_approx=thisdof(k).zbed;
    end
    mea(i)=setfield(mea(i),vname,thismea);
  end
end
clear dof  % no longer needed

% convert sonar data from zobs+data to h
for i=1:length(mea)
  this=mea(i).s;
  for j=1:length(this)
    this(j).data=-(this(j).zobs-this(j).data);
  end
  mea(i).h=this;
end
for i=1:length(mea)
  mea2(i)=rmfield(mea(i),'s');
end
mea=mea2; clear mea2

% extract just the data in the desired time period for this case
for i=1:length(mea)
  if(mea(i).meta.dnum_est<dnum(1))
    istart=i;
  elseif(mea(i).meta.dnum_est>dnum(2))
    iend=i;
    break;
  end
end
mea=mea(istart:iend);

%--------------------------------------
% load forcing conditions: wind, waves
%--------------------------------------

% wave data from nc-file
fn{1}='data/duck94_fulldataset/FRF-ocean_waves_8m-array_199409.nc';
fn{2}='data/duck94_fulldataset/FRF-ocean_waves_8m-array_199410.nc';
waves8m=struct;
waves8m.dnum_est=[];
waves8m.Hrms    =[];
waves8m.sigma   =[];
waves8m.theta0  =[];
waves8m.tide    =[];
for i=1:2
  waves8m.dnum_est=cat(1,waves8m.dnum_est,epoch2Matlab(ncread(fn{i},'time'))-5/24);
  waves8m.Hrms    =cat(1,waves8m.Hrms    ,ncread(fn{i},'waveHs')/1.4);
  waves8m.sigma   =cat(1,waves8m.sigma   ,2*pi./ncread(fn{i},'waveTp'));
  waves8m.theta0  =cat(1,waves8m.theta0  ,deg2rad(71.8-ncread(fn{i},'waveMeanDirection')));
  waves8m.tide    =cat(1,waves8m.tide    ,ncread(fn{i},'waterLevel'));
end

% wind data from nc-file
fn{1}='data/duck94_fulldataset/FRF-met_wind_derived_199409.nc';
fn{2}='data/duck94_fulldataset/FRF-met_wind_derived_199410.nc';
windEOP=struct;
windEOP.dnum_est=[];
windEOP.speed   =[];
windEOP.direc   =[];
for i=1:2
  windEOP.dnum_est=cat(1,windEOP.dnum_est,epoch2Matlab(ncread(fn{i},'time'))-5/24);
  windEOP.speed   =cat(1,windEOP.speed   ,ncread(fn{i},'windSpeed'));
  windEOP.direc   =cat(1,windEOP.direc   ,deg2rad(ncread(fn{i},'windDirection')-71.8));  % cartesian convention
end

% convert wind into vector stress, following Reniers et al. (2004) eqn (6).  Units are N/m2
windEOP.cd=0.002;   % recommended in text
windEOP.rhoa=1;
windEOP.tau(:,1)=windEOP.cd*windEOP.rhoa*windEOP.speed.^2.*cos(windEOP.direc);
windEOP.tau(:,2)=windEOP.cd*windEOP.rhoa*windEOP.speed.^2.*sin(windEOP.direc);

%--------------------------------------
% initial bathymetry: merge SPUV with CRAB
%--------------------------------------

% 1) use 07/08 survey to get the full offshore profile out to x=900m
fn='data/crab/FRF_19940708_0690_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
ind=find(ncread(fn,'profileNumber')==905);
chop=@(data)data(ind);
offshoreprof.x=chop(ncread(fn,'xFRF'));
offshoreprof.y=chop(ncread(fn,'yFRF'));
offshoreprof.z=chop(ncread(fn,'elevation'));

% 2) grab the most-recent survey profile, valid for x<450m or so
data=load(bathyfn);
onshoreprof.x=data(:,1);
onshoreprof.y=data(:,2);
onshoreprof.z=data(:,3);

% 3) grab the initial-time SPUV data for the exact conditions at t=0
this=mea(1).h;
spuvprof.x=[this.x];
spuvprof.y=[this.y];
spuvprof.z=-nanmean([this.data]);

% % v1: merge the data, using loess interp with weights
% x=[offshoreprof.x(:); onshoreprof.x(:); spuvprof.x(:)];
% z=[offshoreprof.z(:); onshoreprof.z(:); spuvprof.z(:)];
% s=[10*ones(size(offshoreprof.z(:)));
%    5*ones(size(onshoreprof.z(:)));
%    1*ones(size(spuvprof.z(:)))];
% lx=100;
% grid.h = -scalecInterp(x,z,s,grid.x,lx,'quadloess');

% v2: merge the data by stitching with a spline
xx=[onshoreprof.x(onshoreprof.x<=min(spuvprof.x)-50);
    spuvprof.x(:);
    offshoreprof.x(offshoreprof.x>=max(spuvprof.x)+50)];
zz=[onshoreprof.z(onshoreprof.x<=min(spuvprof.x)-50);
    spuvprof.z(:);
    offshoreprof.z(offshoreprof.x>=max(spuvprof.x)+50)];
[~,ind]=unique(xx);
grid.h=-interp1(xx(ind),zz(ind),grid.x,'spline');

% do some mild smoothing to remove ripples etc.
m=5;
nx=length(grid.h);
hnew=nan(size(grid.h));
for i=1:nx
  ind=max(1,min(nx,i+[-m:m]));
  hnew(i)=nanmean(grid.h(ind));
end
grid.h=hnew;

% inspect the result
clf, hold on
plot(offshoreprof.x,offshoreprof.z,'k.')
plot(onshoreprof.x,onshoreprof.z,'r.')
plot(spuvprof.x,spuvprof.z,'go')
plot(grid.x,-grid.h,'m-','linewidth',1.5)

% note, hydroSedModel.m requires a grid with +'ve onshore for x.  So, flip
% it around.
grid.h=flipud(grid.h);
grid.xFRF=flipud(grid.x);

%--------------------------------------
% convert observations into the particular format used by assim_1dh.m
%--------------------------------------

% parse the mea struct-array into the variables I care about
for i=1:length(mea)  % loop over 3h sampling intervals
  % disp(['  sampling interval ' num2str(i) ' of ' num2str(length(mea))])
  for n=1:10  % loop over 512s means in 3h interval
    % disp(['    timestamp ' num2str(n) ' of ' num2str(10)])
    thisobs=struct;
    clear dnum_est
    vname='vHh';
    for k=1:length(vname)  % for each data type
      % disp(['      variable ' vname(k)])
      this=getfield(mea(i),vname(k));
      thisdata=struct;
      for j=1:length(this)  % loop over sensor positions
        thisdata.d(j)=this(j).data(n);
        [~,thisdata.ind(j)]=min(abs(grid.xFRF-this(j).x));
        dnum_est{k}(j)=this(j).tstart_est(n);
      end
      thisobs=setfield(thisobs,vname(k),thisdata);
      dnum_est{k}=unique(dnum_est{k});
      if(length(dnum_est{k})>1)
        error('mismatched timestamps')
      end
    end
    thisobs.dnum_est(1)=unique(cell2mat(dnum_est));
    thisobs.v.e=ones(size(thisobs.v.d))*.2;  % v-err
    thisobs.H.e=ones(size(thisobs.H.d))*.1;  % H-err
    thisobs.h.e=ones(size(thisobs.h.d))*.1;  % h-err
    obs(n,i)=thisobs;
  end
end
obs=obs(:);  % time-ordered, 17-26min intervals

%--------------------------------------
% assimilation
%--------------------------------------

verb=1;

% initialize prior sediment transport input parameters and their covariance
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
if(strcmp(sedmodel,'dubarbier'))
  params.Cw=0.00483 ;
  params.Cc=0.02002 ;
  params.Cf=0.01173 ;
  params.Ka=0.631e-4;
  params_covar=diag([0 .001 .005 .005 1e-4].^2);  % TODO, this is a wild guess
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;
  params.m=11;
  params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
  params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
  params_covar=diag([0 .5 2 .1 1].^2);  % TODO, this is a wild guess
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
  params_covar=diag([0 .5 .5].^2);  % TODO, this is a wild guess
else
  error(['invalid sedmodel=' sedmodel])
end

% initialize prior parameter values
clear prior
prior=grid;  % initial grid (x,h,xFRF)
prior.ka_drag=0.0125;  % tuned by Ruessink et al. (2001)
prior.d50=180e-6;
prior.d90=400e-6;
prior.params=params;
prior.sedmodel=sedmodel;

% initialize prior parameter covariances
xx=meshgrid(grid.xFRF);
sig_h=.5*exp(-3*(grid.xFRF-200).^2/150^2);
Ch0=diag(sig_h)*exp(-3*(xx-xx').^2/50^2)*diag(sig_h);
prior.Chs=blkdiag(Ch0,params_covar);
prior.CH0=0.1^2;
prior.Ctheta0=deg2rad(5)^2;
prior.Cka=0.0005^2;

% time loop
for n=1:(length(obs)-1)

  % forcing conditions for this time step
  prior.H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,obs(n).dnum_est);
  prior.theta0=interp1(waves8m.dnum_est,waves8m.theta0,obs(n).dnum_est);
  prior.omega =interp1(waves8m.dnum_est,waves8m.sigma ,obs(n).dnum_est);
  prior.tauw=interp1(windEOP.dnum_est,windEOP.tau,obs(n).dnum_est);
  prior.tauw=repmat(prior.tauw,nx,1);

  % add tide for this time step
  tide=interp1(waves8m.dnum_est,waves8m.tide,obs(n).dnum_est);
  prior.h=prior.h+tide;

  % Define the alongshore pressure gradient detady such that the model will
  % match the measured currents in 8m depth.  Parts of this code are copied
  % from hydro_ruessink2001.m, most notably the equation for the modeled
  % longshore current 'vm'
  rho=1030;
  ind=find(prior.xFRF(obs(n).v.ind)>800);
  if(isempty(ind))
    prior.detady=zeros(nx,1);
  else
    thisx=prior.xFRF(obs(n).v.ind(ind));
    thisv=obs(n).v.d(ind)';
    thish=interp1(prior.xFRF,prior.h,thisx);  % includes tide
    a=1.16;  % empirical constant
    Cd=0.015*(prior.ka_drag./thish).^(1/3);
    g=9.8;
    thisk=zeros(length(thish),1);  % init
    for j=1:length(thish)  % dispersion reln
      thisk(j)=fzero(@(k)prior.omega^2-g*k*tanh(k*thish(j)),.1);
    end
    urms=1.416*prior.H0.*prior.omega./(4*sinh(thisk.*thish));
    vm = @(Fy)sqrt( sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2 ).*sign(-Fy);
    Fy=fminsearch(@(Fy)sum((vm(Fy)-thisv).^2),0);  % choose total forcing to match observed v
    detady = 1/mean(g*thish)*( -Fy + prior.tauw(1,2)/rho );  % solve longshore momentum budget
    prior.detady=detady*ones(nx,1);
  end

  % optional: offshore velocity observations should not be used for assimilation
  ind=find(prior.xFRF(obs(n).v.ind)<800);
  obs(n).v.d  =obs(n).v.d(ind);
  obs(n).v.ind=obs(n).v.ind(ind);
  obs(n).v.e  =obs(n).v.e(ind);

  % Assimilation step: assim_1dh.m Uses the hydrodynamic obs-data to correct
  % the prior bathymetry h at time t(n), then forecasts the bathymetry for
  % the next time step t(n+1)=t(n)+dt.  Remove tide from the results.
  dt=diff([obs(n+[0:1]).dnum_est])*24*60*60;
  posterior=assim_1dh(prior,obs(n),dt,verb);
  posterior.h =posterior.h -tide;
  posterior.hp=posterior.hp-tide;

  % save this posterior for later inspection.  Importantly, it contains
  % predictions of sediment transport parameters for this time step, that
  % aren't carried forward in subsequent timesteps.
  postAll(n)=posterior;

  % Use forecast as prior for the next time step t(n+1)
  prior.h=posterior.hp;
  prior.Chs=posterior.Chsp; % blkdiag(posterior.Chsp(1:nx,1:nx),params_covar);
  prior.params=posterior.params;

  save working.mat
end
