function [obs,bathyobs,grid,waves8m,windEOP]=prepObsData(dnum,bathyfn,duck94Case)
%
% [hydroobs,bathyobs,grid,waves8m,windEOP]=prepObsData(dnum,bathyfn,duck94Case)
%
% Loads data needed for model initialization and hydro-assimilation in duck
% 94 test case.  
%
% INPUTS:
%
% dnum(1),dnum(2) = date limits for model run (EST)
% bathyfn = string pointing to location of CRAB profile data for model initialization
% duck94Case = 'a' 'b' 'c' 'd' or 'e', corresponding to which bar migration
%              "case" is being simulated
%
% OUTPUTS:
%
% hydroobs = array of structs, each element has a datenum (EST), and a set
% of observations that can be assimilated to correct hydrodynamic variables,
% using ./util/assim_1dh.m
%
% bathyobs = struct containing CRAB bathymetry profile data collected during
% the model time period.  This data will be assimilated to adjust sediment
% transport parameters over the entire model run, separate from the
% hydrodynamic data assimilation phase.  The times when the profiles were
% collected are indicated by the 'obsn' field, indeces that point to the
% corresponding time in the hydroobs struct-array.  The reason for doing it
% this way is hydroobs determines the "time steps" for the model run during
% the hydro assimilation phase.
%
% grid = model grid, including initial bathymetry at time dnum(1).  This
% gets defined in the prepObsData function because (a) it uses observational
% data to get the initial bathymetry, including data prepped by this
% function, and (b) the grid is needed in order to define the model indexes
% on which the hydroobs and bathyobs data are defined.
%
% 

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
  if(mea(i).meta.dnum_est<=dnum(1))
    istart=i;
  elseif(mea(i).meta.dnum_est>=dnum(2))
    iend=i;
    break;
  end
end
disp(['starting at ' datestr(mea(istart).meta.dnum_est) ' EST'])
disp(['ending   at ' datestr(mea(iend  ).meta.dnum_est) ' EST'])
mea=mea(istart:iend);

%--------------------------------------
% load forcing conditions: wind, waves
%--------------------------------------

% wave data from nc-file
clear fn
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
clear fn
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

% model grid.  Should not need to change this, except maybe if you want to
% test different grid resolution.
grid.dx=5;
grid.x=[100:grid.dx:900]';
grid.nx=length(grid.x);

% 1) use 07/08 survey to get the full offshore profile out to x=900m
fn='data/crab/FRF_19940708_0690_DUCK94_NAVD88_CRAB_Geodimeter_UTC_v20160513.nc';
ind=find(ncread(fn,'profileNumber')==905);
chop=@(data)data(ind);
offshoreprof.x=chop(ncread(fn,'xFRF'));
offshoreprof.y=chop(ncread(fn,'yFRF'));
offshoreprof.z=chop(ncread(fn,'elevation'));

% 2) For case-a ONLY, grab the initial-time SPUV data for the initial
% bathymetry.  For other cases, use the initial crab profile, filename
% supplied as input parameter 'bathyfn'.  I am using Steve Henderson's
% processed data for SPUV altimeters.
if(duck94Case=='a')
  this=mea(1).h;
  ind=find(925<=[this.y]&[this.y]<=935);  % exclude off-transect sensors
  initdata.x=[this(ind).x];
  initdata.y=[this(ind).y];
  initdata.z=-mean([this(ind).data]);
else
  this=load(bathyfn);
  initdata.x=this(:,1);
  initdata.y=this(:,2);
  initdata.z=this(:,3);
end

% v2: merge the data by stitching with a spline
xx=[offshoreprof.x(offshoreprof.x<=min(initdata.x)-50);
    initdata.x(:);
    offshoreprof.x(offshoreprof.x>=max(initdata.x)+50)];
zz=[offshoreprof.z(offshoreprof.x<=min(initdata.x)-50);
    initdata.z(:);
    offshoreprof.z(offshoreprof.x>=max(initdata.x)+50)];
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

% note, hydroSedModel.m requires a grid with +'ve onshore for x.  So, flip
% it around.
grid.h=flipud(grid.h);
grid.xFRF=flipud(grid.x);

% inspect the result
clf, hold on
plot(offshoreprof.x,offshoreprof.z,'k.')
plot(initdata.x,initdata.z,'go')
plot(grid.xFRF,-grid.h,'m-','linewidth',2)
legend('Offshore section','data t=0','interp result')
% disp('test code'); pause; return;

%--------------------------------------
% convert observations into the particular format used by assim_1dh.m
%--------------------------------------

% parse the mea struct-array into the variables I care about
clear obs
for i=1:length(mea)  % loop over 3h sampling intervals
  % disp(['  sampling interval ' num2str(i) ' of ' num2str(length(mea))])
  for n=1:10  % loop over 512s means in 3h interval
    % disp(['    timestamp ' num2str(n) ' of ' num2str(10)])
    thisobs=struct;
    clear dnum_est
    vname='uvHh';
    for k=1:length(vname)  % for each data type
      % disp(['      variable ' vname(k)])
      this=getfield(mea(i),vname(k));
      thisdata=struct;
      ind=find(925<=[this.y]&[this.y]<=935);  % exclude off-transect sensors
      for j=1:length(ind)  % loop over sensor positions
        if(vname(k)=='h')
          thisdata.d(j)=mean(this(ind(j)).data);  % 3h-average for sonar
        else
          thisdata.d(j)=this(ind(j)).data(n);
        end
        if(vname(k)=='u' | vname(k)=='v')
          thisdata.z(j)=this(ind(j)).zobs-this(ind(j)).zbed_approx;
        end
        [~,thisdata.ind(j)]=min(abs(grid.xFRF-this(ind(j)).x));
        dnum_est{k}(j)=this(ind(j)).tstart_est(n);
      end
      thisobs=setfield(thisobs,vname(k),thisdata);
      dnum_est{k}=unique(dnum_est{k});
      if(length(dnum_est{k})>1)
        error('mismatched timestamps')
      end
    end
    thisobs.dnum_est(1)=unique(cell2mat(dnum_est));
    thisobs.u.e=ones(size(thisobs.u.d))*.2;  % u-err
    thisobs.v.e=ones(size(thisobs.v.d))*.2;  % v-err
    thisobs.H.e=ones(size(thisobs.H.d))*.05;  % H-err
    thisobs.h.e=ones(size(thisobs.h.d))*.05;  % h-err
    obs(n,i)=thisobs;
  end
end
obs=obs(:);  % time-ordered, 17-26min intervals

% remove coincident measurements by averaging
clear obsnew
vname='uvHh';
for i=1:length(obs)
  thisobsnew=struct;
  thisobsnew.dnum_est=obs(i).dnum_est;
  for j=1:length(vname)
    thisdata=getfield(obs(i),vname(j));
    clear thisdatanew
    kk=sort(unique(thisdata.ind));
    for k=1:length(kk)
      thisdatanew.ind(k)=kk(k);
      ind=find(thisdata.ind==kk(k));
      thisdatanew.d(k)=nanmean(thisdata.d(ind));
      thisdatanew.e(k)=nanmean(thisdata.e(ind));
      thisdatanew.ind(k)=nanmean(thisdata.ind(ind));
      if(isfield(thisdatanew,'z'))
        thisdatanew.z(k)=nanmean(thisdata.z(ind));
      end
    end
    thisobsnew=setfield(thisobsnew,vname(j),thisdatanew);
  end
  obsnew(i)=thisobsnew;
end
obs=obsnew;

%--------------------------------------
% Wave asymmetry and skewness observations.  Using SteveH's QC'd SPUV time
% series data.  While I'm at it, also replace u and v with SteveH's data.
% Still use hmo from Elgar's version, since SteveH didn't calculate that but
% Elgar did.
%--------------------------------------

% load spuv data, assign datenums and model indeces
spuv=load('data/SteveH/Duck94all1.mat');
spuv.dnum_est = datenum(1994,9,1,1,0,16):(3/24):datenum(1994,10,19,22,0,16);
for i=1:length(spuv.x)
  if(isnan(spuv.x(i)))
    spuv.ind(i)=nan;
  else
    [dx,ind]=min(abs(spuv.x(i)-grid.xFRF));
    if(dx>5)
      error('should never happen')
    end
    spuv.ind(i)=ind;
  end
end

% get skewness and asymmetry for each u-obs
for i=1:length(obs)
  indn=max(find(spuv.dnum_est<=obs(i).dnum_est));  % 3h chunk that has this time
  dt=(obs(i).dnum_est-spuv.dnum_est(indn))*24*3600;  % time into 3h chunk
  if(dt<0)
    error('I must''ve messed up my timestamping')
  end
  indt=round(dt)*2+[1:2048];  % 2Hz data starting at dt and ending 1024 seconds later
  if(max(indt)>21504)
    error('I didn''t expect an overrun')
    indt=unique(max(21504,indt));  % avoid overrun
  end
  for j=1:length(obs(i).u.ind)
    obs(i).A.ind(j)=obs(i).u.ind(j);
    obs(i).S.ind(j)=obs(i).u.ind(j);
    indx=find(spuv.ind==obs(i).u.ind(j));
    if(length(indx)~=1)
      error(['found ' num2str(length(indx)) ' station matches, should not happen'])
    end
    indt1=setdiff(indt,find(squeeze(isnan(spuv.u(indx,indn,:)))));  % drop nan's
    if(length(indt1)<1024)  % skip station if there are too many nan's
      % disp(['skipping stn ' num2str(j) ' for time ' num2str(i)])
      obs(i).A.d(j)=nan;
      obs(i).A.e(j)=nan;
      obs(i).S.d(j)=nan;
      obs(i).S.e(j)=nan;
      obs(i).u.d(j)=nan;
      obs(i).u.e(j)=nan;
      obs(i).v.d(j)=nan;
      obs(i).v.e(j)=nan;
    else
      uw=detrend(squeeze(spuv.u(indx,indn,indt1)));  % u(t) for this time and station
      dt=0.5;  % sampling rate, seconds
      obs(i).S.d(j)=mean(uw.^3)./mean(uw.^2).^(3/2); % skewness
      uH=imag(hilbert(uw));
      obs(i).A.d(j)=mean(uH.^3)./mean(uH.^2).^(3/2); % asymmetry
      obs(i).u.d(j)=mean(squeeze(spuv.u(indx,indn,indt1)));  % mean x-shore current
      obs(i).v.d(j)=mean(squeeze(spuv.v(indx,indn,indt1)));  % mean longshore current
      obs(i).A.e(j)=0.1;  % error levels are a bit arbitrary
      obs(i).S.e(j)=0.1;  % error levels are a bit arbitrary
      obs(i).u.e(j)=0.1;  % error levels are a bit arbitrary
      obs(i).v.e(j)=0.1;  % error levels are a bit arbitrary
      obs(i).u.z(j)=spuv.h(indx,indn)-spuv.uvpos(indx,indn,3);
    end
  end
  ivalid=find(~isnan(obs(i).S.d));  % drop nan's
  obs(i).S.d=obs(i).S.d(ivalid);
  obs(i).S.e=obs(i).S.e(ivalid);
  obs(i).S.ind=obs(i).S.ind(ivalid);
  ivalid=find(~isnan(obs(i).A.d));
  obs(i).A.d=obs(i).A.d(ivalid);
  obs(i).A.e=obs(i).A.e(ivalid);
  obs(i).A.ind=obs(i).A.ind(ivalid);

  % while I'm at it, replace all h observations with SteveH's version too
  indx=find(~isnan(spuv.h(:,indn)));
  obs(i).h.d=spuv.h(indx,indn);
  obs(i).h.ind=[];
  for j=1:length(indx)
    [~,obs(i).h.ind(j)]=min(abs(spuv.x(indx(j))-grid.xFRF));
  end
  obs(i).h.e=.1*ones(size(obs(i).h.d));

end

clear spuv

%------------------------------------------
% cases b,c,e all have a few bogus Hmo data points involved.  Cull them out.
%------------------------------------------
if(duck94Case=='a' | duck94Case=='d')
  warning('data QC step for H is not tested for cases a or d!!')
end

clear tmp
for i=1:length(obs)
  tmp(i)=max(abs(obs(i).H.d));
end
tmpfilt=medfilt1(tmp);
for i=find(abs(tmpfilt-tmp)>std(tmp)*.5)
  clf, hold on
  plot(obs(i).H.ind,obs(i).H.d,'ko')
  ibad=find(obs(i).H.d > 1.5*mean(obs(i).H.d));
  if(isempty(ibad))
    disp('could not find bad point?  Something is wrong, kicking you into debugger')
    keyboard;
  end
  plot(obs(i).H.ind(ibad),obs(i).H.d(ibad),'ro')
  xlabel('index')
  ylabel('H [m]')
  input(['Removing bad wave height point at time step i=' num2str(i) ', marked in red on plot (press enter)'])
  igod=setdiff(1:length(obs(i).H.ind),ibad);
  obs(i).H.ind=obs(i).H.ind(igod);
  obs(i).H.e=obs(i).H.e(igod);
  obs(i).H.d=obs(i).H.d(igod);
end

%------------------------------------------
% Load data for bathy assimilation phase.  For case-a, there is no CRAB
% survey so SPUV altimeter data are used.  For other cases, load all
% relevant CRAB profiles and decimate.  After this block of code the
% following variables will be defined:
%
% bathyobs(m).h.d =  of observations to be assimilated, for each time m
% bathyobs(m).h.e =  obs-error, for each time m
% measind = list of indeces for obs-data, same for all observation times
% obsn = model time step corresponding to each element in bathyobs
%
%------------------------------------------

if(duck94Case=='a')

  % use SteveH's version of SPUV observations, they are nicely QC'd for
  % bathymetry.  Just grab the observations at the end of the simulation
  % period.
  spuvall=load('data/SteveH/Duck94all1.mat');
  spuvall.dnum_est = datenum(1994,9,1,1,0,0):(3/24):datenum(1994,10,19,22,0,0);
  obsn=length(obs);  % just observe at the last time step
  [dt,n0]=min(abs(spuvall.dnum_est - obs(obsn).dnum_est));
  if(dt>1/24)
    error('can''t find appropriate obs from spuv, should not happen')
  end
  for i=1:length(spuvall.x)
    [~,measind(i)]=min(abs(grid.xFRF-spuvall.x(i)));
  end
  bathyobs.h.d=spuvall.h(:,n0);
  ind=find(~isnan(bathyobs.h.d));
  measind=measind(ind);  % drop nan's
  bathyobs.h.d=bathyobs.h.d(ind);  % drop nan's
  bathyobs.h.e=ones(size(bathyobs.h.d))*.05;
  bathyobs.h.ind=measind;
  bathyobs.dnum_est=obs(obsn).dnum_est;

else

  % CRAB data.  Decimate spatially to keep the number of observations to a
  % reasonable number.  Assume the data are collected at noon each day,
  % coinciding with one of the SPUV collection times (noting the latter
  % determine my model time steps).  Exclude any profiles that don't span
  % a reasonable range in x.
  clear bathyobs obsn
  tmpfn=fileList_ls('data/duck94_fulldataset/Dropbox/Duck94_profiles/*profile');
  tmpx=650:5:800;  % define obs points in grid coords (m)
  [~,measind]=intersect(grid.x',tmpx);
  for i=1:length(tmpfn)
    tmpdata=load(tmpfn{i});
    this=struct;
    this.h.ind=measind;
    [~,ii]=unique(tmpdata(:,1));
    this.h.d=interp1(tmpdata(ii,1),-tmpdata(ii,3),grid.xFRF(measind));
    this.h.e=.05*ones(size(measind));  % obs error
    tmpdstr=strsh(tmpfn{i},'t');
    tmpdnum=datenum(1994,str2num(tmpdstr(1:2)),str2num(tmpdstr(3:4)),12,0,0);
    [dt,ind]=min(abs([obs.dnum_est]-tmpdnum));
    if(abs(dt)<.1 & sum(isnan(this.h.d))==0)
      this.dnum_est=obs(ind).dnum_est;
      % disp(['adding crab survey ' datestr(this.dnum_est)])
      if(~exist('bathyobs'))  % first one
        bathyobs=this;
        obsn=ind;
      else  % append to list
      bathyobs(end+1)=this;
      obsn(end+1)=ind;
      end
    else  % survey not in range of model simulation
      ; % disp(['SKIPPING crab survey ' datestr(tmpdnum) ' (dt=' num2str(dt) ')'])
    end
  end
  clear tmp* this

end

% wrap everything up for output
for n=1:length(bathyobs)
  bathyobs(n).obsn=obsn(n);
  bathyobs(n).measind=measind;
end

% Note, for everything except case-e, will just assimilate the "final"
% bathymetry to save time, drop the rest of them
if(duck94Case~='e')
  disp('only keeping the final bathy profile')
  bathyobs=bathyobs(end);
else
  disp(['case-e: Will assimilate ' num2str(length(bathyobs)) ' available bathymetry profiles, only dropping the initial bathymetry.'])
  bathyobs=bathyobs(2:end);
end

% Optional check: Do the bathy data look reasonable??  It should show a nice
% clean bar migration, else I need to curate the data a bit to ensure I only
% assimilate the good stuff.
clf, hold on
lstr={};
cc=interp1(1:64,jet(64),linspace(1,64,length(bathyobs)));
for n=1:length(bathyobs)
  plot(grid.xFRF(bathyobs(n).h.ind),-bathyobs(n).h.d,'o-','color',cc(n,:))
  lstr{n}=datestr(bathyobs(n).dnum_est);
end
legend(lstr)
pause(.1);
% input(['Visual check if bathy data are clean enough for assimilation.  Press enter to continue.']);
