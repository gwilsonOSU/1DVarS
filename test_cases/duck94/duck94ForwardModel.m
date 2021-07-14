% The original version of this code is saved in the same directory. This is
% just a branch for experiments and debugging.

%NOTE: The change log also keeps track of the changes made to all the
%functions called by this script for convenience.

% Test cases of forward model skill during canonical on/offshore migration
% cases from Duck94.
%
% The code is largely based on a previous version:
% ./initialExperiments/duck94Assimilation.m.  Some functionality for
% assimilation is included here, however.  See below.
%
% BASIC PURPOSE: Forward modeling of selected time periods in Duck94 that
% had significant bar migration (see user inputs to select a time period).
% You can choose from one of Dubarbier, van der A, or Soulsby Van-Rijn
% sediment transport formulations, and compare their skill for predicting
% the bar migration.
%
% OPTIONS / EXPERIMENTS:
%
% See user inputs block at top of code, you can turn on/off the following:
%
% a) Assimilate data to correct hydrodynamics (waves and currents).
%    Doing so results in a much better sediment transport prediction.
%
% b) Include adjoint sensitivity analysis after each time step.  This could
%    be used to look for the most-sensitive model input parameters w.r.t.
%    morphology change.
%
% addpath(genpath('../../src'))  % hydroSedModel.m, and dependencies
% addpath ./data/tools
clear

%--------------------------------------
% user inputs
%--------------------------------------

% optional flags
doassim=1;   % use assimilation to correct hydro errors
doadjoint=0; % include a sensitivity analysis for sed transport inputs

% duck94 case: set this to a,b,c, or d...  These follow Gallagher et
% al. (1998) four standard test cases.  All but case (b) are offshore bar
% migration events.  Note, case-b is a good example of an onshore migration,
% and currently (last checked June, 2021) this case gives good results with
% options sedmodel='vanderA' and doassim=1.
duck94Case='b';


%Adding directory to save outputs:
dir = strcat('/home/');
if exist ('dir','var') == 0
    dir = pwd;
    mkdir Ch
    mkdir ChPlots
    mkdir AssimilatedBathymetry
end





% sediment model: uncomment one of the following...
% sedmodel='dubarbier';
sedmodel='vanderA';

% model grid.  Should not need to change this, except maybe if you want to
% test different grid resolution.
grid.dx=5;
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
  warning('test code, changed start time'); dnum(1)=datenum(1994,9,24,19,0,0);  % TEST CODE, start at a time with a nice bar profile
  % bathyfn='data/crab/FRF_19940921_0710_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0921_profile';
 case 'c'
  dnum(1)=datenum(1994,10,2,22,0,0);  % Gallagher et al. used 1900-2200, but
                                      % no shoreface sensors available then?
  dnum(2)=datenum(1994,10,4,16,0,0);
  warning('test code, changed start time'); dnum(1)=728570.3;  % TEST CODE, start at a time with non-trivial waves
  % bathyfn='data/crab/FRF_19941004_0716_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1003_profile';
 case 'd'
  dnum(1)=datenum(1994,10,10,22,0,0);
  dnum(2)=datenum(1994,10,15,22,0,0);
  % bathyfn='data/crab/FRF_19941010_0722_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1021_profile';
end

%--------------------------------------
% Get SPUV measurements for this time period.  Note, code to read the raw
% data is included but is commented out.  Instead, I saved a local mat-file
% version of the data to speed up the loading process.
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
waves8m.sigmap  =[];
waves8m.sigmam  =[];
waves8m.theta0  =[];
waves8m.tide    =[];
for i=1:2
  waves8m.dnum_est=cat(1,waves8m.dnum_est,epoch2Matlab(ncread(fn{i},'time'))-5/24);
  waves8m.Hrms    =cat(1,waves8m.Hrms    ,ncread(fn{i},'waveHs')/1.4);
  waves8m.sigmap  =cat(1,waves8m.sigmap  ,2*pi./ncread(fn{i},'waveTp'));
  waves8m.sigmam  =cat(1,waves8m.sigmam  ,2*pi./ncread(fn{i},'waveTm1'));
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
nx=length(grid.h);

% % v1: MA-smoothing to remove ripples etc.
% m=5;
% hnew=nan(size(grid.h));
% for i=1:nx
%   ind=max(1,min(nx,i+[-m:m]));
%   hnew(i)=nanmean(grid.h(ind));
% end
% grid.h=hnew;

% v2: low-pass filter to remove ripples, etc.
fc = 1/50;   % band lower bound in 1/m
fs=1/abs(grid.x(2)-grid.x(1));
[b,a] = butter(5,fc/(fs/2),'low');
htmp = detrend(grid.h);
htrend = grid.h-htmp;
htmpA = filtfilt(b,a,htmp);
hfilt = htmpA + htrend;
grid.h=hfilt;

% % inspect the result
% clf, hold on
% plot(offshoreprof.x,offshoreprof.z,'k.')
% plot(onshoreprof.x,onshoreprof.z,'r.')
% plot(spuvprof.x,spuvprof.z,'go')
% plot(grid.x,-grid.h,'m-','linewidth',1.5)

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
% Model initialization
%--------------------------------------

% initialize sediment transport input parameters
params=struct;
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
if(strcmp(sedmodel,'dubarbier'))
  params.Cw=.00483;  % default 0.000483;  % hsu et al. 0.0046
  params.Cc=.02002 ;  % default 0.02002;  % hsu et al. 0.0053
  params.Cf=0.01173;  % default 0.01173 ; hsu et al. 0.0
  params.Ka=0.631e-4;  % default 0.631e-4;  % hoefel & elgar 1.4e-4; hsu et al. 0.0
%   params_std=[.0002 .01  .005 5e-5];   % Old line
  params_std=[5e-3 .0002 .01  .005 5e-5];
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;  % vanderA 1.2.  Larger values promote offshore bar migration
  params.m=11;  % vanderA 11; hsu et al. 11.  Just scales everything up
  params.xi=1.7;  % ??? tuning parameter, O(1) according to Kranenburg (2013).  VDA pg. 35 says 1.7
  params.alpha=20;  % comes in eqn 27-28, not the same as eqn 19
  params_std=[5e-3 .2 2 .5 2];
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% initialize 'modelinput' struct.  This will be used as the base input for
% each model simulation.  Some of the inputs will get updated with each time
% step, others will remain constant for all time steps.
modelinput=grid;  % initialize: x,h,xFRF
modelinput.ka_drag=0.0125;  % tuned by Ruessink et al. (2001)
modelinput.d50=nan(nx,1);
modelinput.d50(grid.xFRF<150)=400e-6;  % shoreface coarse sand, Birkemeier et al. (1985)
modelinput.d90(grid.xFRF<150)=400e-6;  % shoreface coarse sand
modelinput.d50(grid.xFRF>=150)=180e-6; % Dubarbier et al. (2015) used 170um here, but 180 is from Birkemeier et al. (1985)
modelinput.d90(grid.xFRF>=150)=240e-6;  % Birkemeier et al., 1985
modelinput.params=params;
modelinput.sedmodel=sedmodel;

% initialize covariances for assimilation step
xx=meshgrid(grid.xFRF);
sig_h=0.01*exp(-3*(grid.xFRF-200).^2/150^2);  % minimal h-correction
Ch0=diag(sig_h)*exp(-3*(xx-xx').^2/50^2)*diag(sig_h);
modelinput.Cgamma=.1^2*exp(-3*(xx-xx').^2/25^2);
modelinput.Chs=blkdiag(Ch0,diag(params_std.^2));
modelinput.CH0=0.25^2;
modelinput.Ctheta0=deg2rad(10)^2;
modelinput.Cka=0.005^2;

%--------------------------------------
% Time loop of model runs
%--------------------------------------

%Declaring the Variable to store Ch values:
InputChStore = nan (length(obs),1);
%Declaring the Variable to store rmse values:
RMSE = nan (length(obs),1);
RMSEp = nan (length(obs),1);

%Declaring the time variable:
Time = zeros (length(obs),1);



for n=1:(length(obs)-1)  % note: recommend to start at n=140 for case-b testing
  disp(['time ' num2str(n) ' of ' num2str(length(obs))])

  % Set wave conditions for this time step.  Note, use Tm1 wave period
  % (corresponding to waves8m.sigmam) following Ruessink et al. (2012)
  % formulation for wave shape
  modelinput.H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,obs(n).dnum_est);
  modelinput.theta0=interp1(waves8m.dnum_est,waves8m.theta0,obs(n).dnum_est);
  modelinput.omega =interp1(waves8m.dnum_est,waves8m.sigmam,obs(n).dnum_est);
  modelinput.tauw=interp1(windEOP.dnum_est,windEOP.tau,obs(n).dnum_est);
  modelinput.tauw=repmat(modelinput.tauw,nx,1);
  modelinput.dgamma=zeros(nx,1);

  % Add tide for this time step.  Note by convention modelinput.h does not
  % contain tide.  We add tide before each forecast step (now) and remove it
  % after the forecast is done.
  tide=interp1(waves8m.dnum_est,waves8m.tide,obs(n).dnum_est);
  modelinput.h=modelinput.h+tide;

  % define observations for assimilation
  thisobs=obs(n);
  thisobs.h.d=thisobs.h.d+tide;
  thisobs=rmfield(thisobs,'v');  % don't assimilate v
  thisobs=rmfield(thisobs,'h');  % don't assimilate h
  thisobs.H.e=.05*ones(size(thisobs.H.e));

  % Define the alongshore pressure gradient detady such that the model will
  % match the measured currents in 8m depth.  Parts of this code are copied
  % from hydro_ruessink2001.m, most notably the equation for the modeled
  % longshore current 'vm'
  rho=1030;
  ind=find(modelinput.xFRF(obs(n).v.ind)>800);
  if(isempty(ind))
    modelinput.detady=zeros(nx,1);
  else
    thisx=modelinput.xFRF(obs(n).v.ind(ind));
    thisv=obs(n).v.d(ind)';
    thish=interp1(modelinput.xFRF,modelinput.h,thisx);  % includes tide
    a=1.16;  % empirical constant
    Cd=0.015*(modelinput.ka_drag./thish).^(1/3);
    g=9.8;
    thisk=zeros(length(thish),1);  % init
    for j=1:length(thish)  % dispersion reln
      thisk(j)=fzero(@(k)modelinput.omega^2-g*k*tanh(k*thish(j)),.1);
    end
    urms=1.416*modelinput.H0.*modelinput.omega./(4*sinh(thisk.*thish));
    vm = @(Fy)sqrt( sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2 ).*sign(-Fy);
    Fy=fminsearch(@(Fy)sum((vm(Fy)-thisv).^2),0);  % choose total forcing to match observed v
    detady = 1/mean(g*thish)*( Fy - modelinput.tauw(1,2)/rho );  % solve longshore momentum budget
    modelinput.detady=detady*ones(nx,1);
  end

  % Run model to step forward in time.  This takes the initial 'modelinput'
  % struct, and produces a new 'fcst' struct containing a forecast valid for
  % the current time step.
  dt=diff([obs(n+[0:1]).dnum_est])*24*60*60;
  mfact=1;
  dtm=dt/mfact;  % sub-divide into mfact time steps to avoid instability
  
  %Saving the input Ch (sum of its diagonal members)
  InputChStore(n,1) = trace(modelinput.Chs);
  
   % Saving Ch Matrices and Plots
    ChFullMatrix = modelinput.Chs(1:nx, 1:nx); 
    save (strcat(dir,'Ch/','ChFull',num2str(n),'.mat'),'ChFullMatrix')
    
    figure('Units', 'inches','Position', [0,0,8,4])
    s = pcolor(grid.xFRF,grid.xFRF,ChFullMatrix);
    title(strcat('Ch for n = ',num2str(n)),'FontSize',14,'Interpreter','Latex')
    s.FaceColor = 'interp';
    shading flat
    c = colorbar;
    dc = (max(max(ChFullMatrix)) - min(min(ChFullMatrix)))/10;
    c.Ticks = [min(min(ChFullMatrix)): dc: max(max(ChFullMatrix))];
    colormap jet     
    print ('-dpng', '-r300', strcat(dir,'ChPlots/','Ch- ','n=', num2str(n), ' of ', num2str(length(obs)),'.png'))
    close
  
  
  for nn=1:mfact
    if(mfact>1)
      disp(['  sub-step ' num2str(nn) ' of ' num2str(mfact)])
    end

    % Forecast step.  Depending on user input flag 'doassim', this can either be
    % a pure forecast for t=t(n)+dtm, or it can be an assimilative update
    % for t(n) followed by a forecast for t(n)+dtm.  For the latter case,
    % see README for a description of the forecast-assimilation scheme.
    try
    if(~doassim) % v1: without assimilation
      fcst=hydroSedModel(modelinput.x,modelinput.h,modelinput.H0,modelinput.theta0,modelinput.omega,modelinput.ka_drag,modelinput.tauw,modelinput.detady,modelinput.dgamma,modelinput.d50,modelinput.d90,modelinput.params,sedmodel,dtm);
    else % v2: assimilate data to correct hydro model
      [posterior,fcst]=assim_1dh2(modelinput,thisobs,dt,0,1);          
    end
    catch
        warning(strcat('Model crashed at n = ',num2str(n)));
        ErrorMsg = strcat('It crashed at n = ',num2str(n));
        save (strcat(dir,'ErrorMsg.mat'),'ErrorMsg')
        exit
    end


    modelinput.h=fcst.hp;
  end

  % Remove tide again.  We want to keep h in navd88 except when running the
  % forecast step.
  fcst.h =fcst.h -tide;
  fcst.hp=fcst.hp-tide;

  % save the result in an array
  out(n)=fcst;

  % plot results for this time step
  lw=1.5;
  clf
  subplot(221), hold on
  set(gcf,'visible','off')
  plot(grid.xFRF,grid.h,'g','linewidth',lw)
  plot(grid.xFRF,out(n).h,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).h.ind),obs(n).h.d,'ko')
  plot(grid.xFRF,modelinput.h-tide -sqrt(diag(modelinput.Chs(1:nx,1:nx))),'r--')
  plot(grid.xFRF,modelinput.h - tide +sqrt(diag(modelinput.Chs(1:nx,1:nx))),'r--')
  set(gca,'ydir','r')
  legend('initial, n=1',['n=' num2str(n) ' of ' num2str(length(obs))],'observed')
  ylabel('h [m]')
  subplot(222), hold on
  plot(grid.xFRF,out(n).Hrms,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).H.ind),obs(n).H.d,'ko')
  ylabel('H_{rms} [m]')
  subplot(223), hold on
  plot(grid.xFRF,out(n).Q,'b','linewidth',lw)  % out(n).hp-modelinput.h+tide
  ylabel('Q [m^2/s]')
  ylim([-1 1]*5e-4)
  subplot(224), hold on
  plot(grid.xFRF,out(n).vbar,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).v.ind),obs(n).v.d,'ko')
  ylabel('v [m/s]')
  for j=1:4
    subplot(2,2,j)
    box on
    xlim([100 400])
    ax=axis;
    if(j==1)
      [~,ishore]=min(abs(out(n).h+tide));
      plot(ax(1:2),-[1 1]*tide,'k--')
    else
      plot(ax(1:2),[0 0],'k--')
    end
    plot(grid.xFRF(ishore)*[1 1],ax(3:4),'k--')
    axis(ax)
%     grid on           %grid is declared as a struct earlier in the script
  end
  subplot(221)
  title(['yday ' num2str(obs(n).dnum_est-datenum(1994,1,0)) ...
        ', ' datestr(obs(n).dnum_est) 'EST'])
  %Saving the plot:
  print ('-dpng', '-r300', strcat(dir,'AssimilatedBathymetry/','Bathymetry- ','n=', num2str(n), ' of ', num2str(length(obs)),'.png'))    
  pause(.01)

  % Optional: having done the forward model simulation, now do the adjoint
  % sensitivity analysis for the sediment transport model.
  if(doadjoint)
    if(~strcmp(sedmodel,'vanderA'))
      error('need to implement adjoint calculation for models other than vanderA')
    end
    eparam=1;  % set to 1 to get sensitivity for each gridpoint, not summed over gridpoints
    [ad_d50(:,n),ad_d90(:,n),ad_h(:,n),ad_Hrms(:,n),ad_kabs(:,n),ad_omega(n),ad_udelta(:,:,n),ad_ws(:,n),ad_param(:,n)] = ...
     ad_qtrans_vanderA(ones(nx,1),fcst.bkgd_qtrans,eparam);
  end


% %  %CALCULATING rmse
    RMSE(n) = rms(out(n).h(obs(n).h.ind) - obs(n).h.d');
    RMSEp(n+1) = rms(out(n).hp(obs(n+1).h.ind) - obs(n+1).h.d');
    
% Storing time:
    if n < length(obs)      %Using the if condition to avoid crashing if the outer loop is expanded to length(obs) instead of length(obs) -1
        Time(n+1) = Time(n) + hours(datetime(obs(n+1).dnum_est,'ConvertFrom','datenum') - datetime(obs(n).dnum_est,'ConvertFrom','datenum'));
    end
  
  
  
  % Use the forecast bathymetry as input for the next time step t(n+1)
  modelinput.h=out(n).hp;
  modelinput.Chs(1:nx, 1:nx) = posterior.Chsp(1:nx, 1:nx);
  
  %Saving Ch, RMSE, Time, Posterior
save (strcat(dir,'InputChLog.mat'),'InputChStore')
save (strcat(dir,'RMSE.mat'),'RMSE')
save (strcat(dir,'RMSEp.mat'),'RMSEp')
save (strcat(dir,'Time.mat'),'Time')

SaveOut = out(n);
save (strcat(dir,'Output/out', num2str(n),'.mat'),'SaveOut')
save (strcat(dir,'Posterior/posterior', num2str(n),'.mat'),'posterior')
end



