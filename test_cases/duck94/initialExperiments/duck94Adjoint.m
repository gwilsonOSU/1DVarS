% Test cases of adjoint sensitivity during canonical on/offshore migration
% cases from Duck94.
%
% This is mostly copied from duck94Forwardmodel.m
%
addpath(genpath('~/work/unfunded_projects/sedimentTransport1D_TLAD/hydroSedModel/'))
addpath tools
clear

%--------------------------------------
% user inputs
%--------------------------------------

% duck94 case: set this to a,b,c, or d...  These follow Gallagher et
% al. (1998) four standard test cases.  All but case (b) are offshore bar
% migration events.  Set 'casen' to the time step to be simulated with the
% adjoint model.
duck94Case='b';
casen=200;  % time step at which to do adjoint sensitivty analysis

% sediment model: uncomment one of the following... note, for adjoint
% modeling I've only coded the vanderA case, not that it would be hard to do
% the others
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';
if(~strcmp(sedmodel,'vanderA'))
  error('models other than vanderA not supported in this version, but could be added easily')
end

% model grid
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
% model run
%--------------------------------------

% initialize sediment transport input parameters
params=struct;
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
if(strcmp(sedmodel,'dubarbier'))
  params.Cw=.046;  % default 0.000483;  % hsu et al. 0.0046
  params.Cc=.053 ;  % default 0.02002;  % hsu et al. 0.0053
  params.Cf=0;  % default 0.01173 ; hsu et al. 0.0
  params.Ka=0;  % default 0.631e-4;  % hoefel & elgar 1.4e-4; hsu et al. 0.0
  params_std=[.0002 .01  .005 5e-5];
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;  % vanderA 1.2.  Larger values promote offshore bar migration
  params.m=11;  % vanderA 11; hsu et al. 11.  Just scales everything up
  params.xi=1.7;  % ??? tuning parameter, O(1) according to Kranenburg (2013).  VDA pg. 35 says 1.7
  params.alpha=8.2;  % comes in eqn 27-28, not the same as eqn 19
  params_std=[.2 2 .5 2];
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% initialize 'model' struct
clear model
model=grid;  % initial grid (x,h,xFRF)
model.ka_drag=0.0125;  % tuned by Ruessink et al. (2001)
model.d50=nan(nx,1);
model.d50(grid.xFRF<150)=400e-6;
model.d50(grid.xFRF>=150)=180e-6;  % Dubarbier et al., 2015 used 170um here
model.d90=400e-6;
model.params=params;
model.sedmodel=sedmodel;

% initialize covariances for assimilation step
xx=meshgrid(grid.xFRF);
sig_h=0.01*exp(-3*(grid.xFRF-200).^2/150^2);  % minimal h-correction
Ch0=diag(sig_h)*exp(-3*(xx-xx').^2/50^2)*diag(sig_h);
model.Cgamma=.1^2*exp(-3*(xx-xx').^2/25^2);
model.Chs=blkdiag(Ch0,diag(params_std.^2));
model.CH0=0.25^2;
model.Ctheta0=deg2rad(10)^2;
model.Cka=0.005^2;

% rather than a time loop, here I just select a particular time in the simulation
n=casen;
  disp(['time ' num2str(n) ' of ' num2str(length(obs))])

  % forcing conditions for this time step.  Note, use Tm1 period following
  % Ruessink et al. (2012) formulation for wave shape
  model.H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,obs(n).dnum_est);
  model.theta0=interp1(waves8m.dnum_est,waves8m.theta0,obs(n).dnum_est);
  model.omega =interp1(waves8m.dnum_est,waves8m.sigmam,obs(n).dnum_est);
  model.tauw=interp1(windEOP.dnum_est,windEOP.tau,obs(n).dnum_est);
  model.tauw=repmat(model.tauw,nx,1);
  model.dgamma=zeros(nx,1);

  % add tide for this time step
  tide=interp1(waves8m.dnum_est,waves8m.tide,obs(n).dnum_est);
  model.h=model.h+tide;

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
  ind=find(model.xFRF(obs(n).v.ind)>800);
  if(isempty(ind))
    model.detady=zeros(nx,1);
  else
    thisx=model.xFRF(obs(n).v.ind(ind));
    thisv=obs(n).v.d(ind)';
    thish=interp1(model.xFRF,model.h,thisx);  % includes tide
    a=1.16;  % empirical constant
    Cd=0.015*(model.ka_drag./thish).^(1/3);
    g=9.8;
    thisk=zeros(length(thish),1);  % init
    for j=1:length(thish)  % dispersion reln
      thisk(j)=fzero(@(k)model.omega^2-g*k*tanh(k*thish(j)),.1);
    end
    urms=1.416*model.H0.*model.omega./(4*sinh(thisk.*thish));
    vm = @(Fy)sqrt( sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2 ).*sign(-Fy);
    Fy=fminsearch(@(Fy)sum((vm(Fy)-thisv).^2),0);  % choose total forcing to match observed v
    detady = 1/mean(g*thish)*( Fy - model.tauw(1,2)/rho );  % solve longshore momentum budget
    model.detady=detady*ones(nx,1);
  end

  % run model to step forward in time
  dt=diff([obs(n+[0:1]).dnum_est])*24*60*60;
  mfact=1;
  dtm=dt/mfact;  % sub-divide into mfact time steps to avoid instability
  for nn=1:mfact
    if(mfact>1)
      disp(['  sub-step ' num2str(nn) ' of ' num2str(mfact)])
    end

    % % v1: without assimilation
    % fcst=hydroSedModel(model.x,model.h,model.H0,model.theta0,model.omega,model.ka_drag,model.tauw,model.detady,model.dgamma,model.d50,model.d90,model.params,sedmodel,dtm);

    % v2: assimilate data to correct hydro model
    [~,fcst]=assim_1dh(model,thisobs,dt,0,0);

    model.h=fcst.hp;
  end

  % don't put the tide back in, because I want to do adjoint modeling next
  % fcst.h =fcst.h -tide;
  % fcst.hp=fcst.hp-tide;

  % show results
  lw=1.5;
  clf
  subplot(221), hold on
  plot(grid.xFRF,grid.h+tide,'g','linewidth',lw)  % note, add tide b/c it has not been removed from fcst
  plot(grid.xFRF,fcst.hp,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).h.ind),obs(n).h.d,'ko')
  set(gca,'ydir','r')
  legend('initial, n=1',['n=' num2str(n) ' of ' num2str(length(obs))],'observed')
  ylabel('h [m]')
  subplot(222), hold on
  plot(grid.xFRF,fcst.Hrms,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).H.ind),obs(n).H.d,'ko')
  ylabel('H_{rms} [m]')
  subplot(223), hold on
  plot(grid.xFRF,fcst.Q,'b','linewidth',lw)  % fcst.hp-model.h+tide
  ylabel('Q [m^2/s]')
  ylim([-1 1]*1e-4)
  subplot(224), hold on
  plot(grid.xFRF,fcst.vbar,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).v.ind),obs(n).v.d,'ko')
  ylabel('v [m/s]')
  for j=1:4
    subplot(2,2,j)
    box on
    xlim([100 400])
    ax=axis;
    if(j==1)
      [~,ishore]=min(abs(fcst.h+tide));
      plot(ax(1:2),-[1 1]*tide,'k--')
    else
      plot(ax(1:2),[0 0],'k--')
    end
    plot(grid.xFRF(ishore)*[1 1],ax(3:4),'k--')
    axis(ax)
    grid on
  end
  subplot(221)
  title(['yday ' num2str(obs(n).dnum_est-datenum(1994,1,0)) ...
        ', ' datestr(obs(n).dnum_est) 'EST'])
  pause(.01)

% having done the forward model simulation, now do the adjoint sensitivity
% analysis for the transport model.  TODO, add other models besides vanderA
% as options.  Note, use the 'eparam' option here to get parameter
% sensitivity at each gridpoint instead of summed over gridpoints.
[ad_d50,ad_d90,ad_h,ad_Hrms,ad_kabs,ad_omega,ad_udelta,ad_ws,ad_param] = ...
    ad_qtrans_vanderA(ones(nx,1),fcst.bkgd_qtrans,1);

% check magnitude of terms for a specific gridpoint
i=133;
disp([' ad_d50(i) = '         num2str(ad_d50   (i)*fcst.d50   (i))]);
disp([' ad_d90(i) = '         num2str(ad_d90*fcst.d90   )]);
disp([' ad_h(i) = '           num2str(ad_h     (i)*fcst.h     (i))]);
disp([' ad_Hrms(i) = '        num2str(ad_Hrms  (i)*fcst.Hrms  (i))]);
disp([' ad_kabs(i) = '        num2str(ad_kabs  (i)*fcst.kabs  (i))]);
disp([' ad_omega(i) = '       num2str(ad_omega*fcst.omega )]);
disp([' ad_udelta(i) = '      num2str(ad_udelta(i)*fcst.udelta(i))]);
disp([' ad_ws(i) = '          num2str(ad_ws    (i)*fcst.ws    (i))]);
disp([' ad_param(i).n = '     num2str(ad_param (i).n    *fcst.params.n    )])
disp([' ad_param(i).m = '     num2str(ad_param (i).m    *fcst.params.m    )])
disp([' ad_param(i).xi = '    num2str(ad_param (i).xi   *fcst.params.xi   )])
disp([' ad_param(i).alpha = ' num2str(ad_param (i).alpha*fcst.params.alpha)])

% plot the result, normalizing adjoints by their bkgd values.  Exclude ws,
% xi, and alpha, because in this case the phase lag effect is not active
scalefact=.1;
clf
subplot(311)
plot(grid.xFRF,fcst.h)
ylabel('h [m]')
set(gca,'ydir','r')
subplot(312)
plot(grid.xFRF,fcst.Q)
ylabel('Q [m^2/s]')
vname={'d50','d90','h','Hrms','kabs','omega','udelta'};%,'ws'};
vname_param={'n','m'};%,'xi','alpha'};
subplot(313), hold on
for i=1:length(vname)
  plot(grid.xFRF,eval(['ad_' vname{i} '.*fcst.' vname{i} '*scalefact']))
end
for i=1:length(vname_param)
  plot(grid.xFRF,eval(['[ad_param.' vname_param{i} '].*fcst.params.' vname_param{i} '*scalefact']))
end
ylabel('dQ [m^2/s] per 10% increase in parameter')
legend(cat(2,vname,vname_param))
ylim([-1 1]*5e-5)
