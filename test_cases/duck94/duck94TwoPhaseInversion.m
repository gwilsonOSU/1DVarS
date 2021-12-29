% Main function for Duck94 assimilative tests, for the four bar migration
% events in the Duck94, which are referred to here as case-a, case-b,
% case-c, and case-d.  All but case-b are offshore bar migration events.
%
% The approach used by this code is a "two-phase" assimilation.  In phase-1,
% hydrodynamic model outputs are refined by assimilating hydrodynamic
% observations.  The result of phase-1 is a background model simulation that
% has very good hydrodynamics, along with a corresponding prediction of
% bathymetry evolution over time.  In phase-2, a measurement of the
% bathymetry at the end of the model period is used to correct sediment
% transport parameters, in order to improve the morphodynamic prediction.
% The result of phase-2 is a set of updated sediment transport parameters
% that are expected to produce a better model outcome.  Phase-1 is then
% repeated, using the new parameters that were just found in phase-2, and
% making any further needed hydrodynamic adjustments.  The two-phase process
% is iterated repeatedly until convergence.
%
% The inputs to this code include the Duck94 bar migration event you want to
% simulate (can be a,b,c,d), and the sediment transport formulation you want
% to use (either van Der A et al. (2013), or Dubarbier et al. which is based
% on the extended energetics model of Hsu et al. (2006)).
%
% This main function uses several sub-functions that help to modularize the
% code and hopefully make things easier to follow.  It generates some
% output .mat files of the results.
%
addpath(genpath('../../../src'))  % hydroSedModel.m and its dependencies
addpath util  % helper functions
% clear

%--------------------------------------
% user inputs
%--------------------------------------

% duck94 case: set this to a,b,c, or d
% duck94Case='d';

% sediment model: uncomment one of the following...
% sedmodel='dubarbier';
sedmodel='vanderA';

% number of two-phase assimilation iterations
nitermax=6;

% output options
dosave=1;  % save figures and .mat files for each iteration

%--------------------------------------
% End user inputs.  Do not edit beyond this point if you are just running cases.
%--------------------------------------

% Set case-dependent model-run parameters, depending on user input
% 'duck94Case'
%
%   dnum(1) = EST datenum for start of this case
%   dnum(2) = EST datenum for end of this case
%   bathyfn = initial bathymetry profile data set
%
% Most of the start/end times are based on Rafati et al. (2021), who
% sensibly chose time periods based on existence of CRAB initial and final
% profile data.  Any deviations from Rafati et al. (2001) are explained in
% the comments.
%
% Note the Duck94 CRAB surveys were always at about noon EST, so I will use
% that as my start/end time-of-day in all cases.
%
switch duck94Case
 
 case 'a'
 
  % NOTE: case-a is special case, cannot use CRAB bathy for initialization!
  % Hence there is no bathyfn defined here.  The initial bathy will instead
  % be handled as a special case later on in this code.
  dnum(1)=datenum(1994,9,1,12,0,0);
  dnum(2)=datenum(1994,9,5,12,0,0);
  bathyfn='';

 case 'b'

  % % Rafati et al. started case-b on the 21st, ended on the 26th.  BUT, see
  % % below notes/tweaks.
  % dnum(1)=datenum(1994,9,21,12,0,0);
  % dnum(2)=datenum(1994,9,26,12,0,0);
  % bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0921_profile';

  % Dubarbier et al. skipped past the wave event on the 21st, and instead
  % defined case-b as starting on Sep 24 when a crab profile was taken.
  % Similarly, Hsu et al. (2006) started their case on the 22nd (though no
  % crab profile is available so I have to start later)
  dnum(1)=datenum(1994,9,24,12,0,0);
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0924_profile';

  % Rafati et al. ended their case-b on sep26, which seems odd.  There is
  % another full crab profile on sep30 that shows further onshore bar
  % migration, and ought to be included in the simulation.  Note Hsu et
  % al. (2006) used Sep27 as their end time for this case, though there is
  % no CRAB survey on that day.  For the above reasons, I will use sep30 as
  % the end date.
  dnum(2)=datenum(1994,9,30,12,0,0);

 case 'c'

  dnum(1)=datenum(1994,9,30,12,0,0);
  dnum(2)=datenum(1994,10,4,12,0,0);
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0930_profile';

  % OPTIONAL: case-c has a long quiescent period at the beginning.  For
  % testing purposes, just skip to the storm period.
  warning('saving time by starting case-c at oct3,0000EST. Be sure the obsdataCache reflects this!');
  dnum(1)=datenum(1994,10,3,0,0,0);

 case 'd'

  dnum(1)=datenum(1994,10,10,12,0,0);
  dnum(2)=datenum(1994,10,14,12,0,0);
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1010_profile';

end

% Load the model grid and grid-referenced observational data for this case.
% Cached as mat-files to save time.
obsdatafn=['obsdataCache/obsdata_case' duck94Case '.mat'];
if(~isempty(dir(obsdatafn)))
  disp(['loading pre-cached obsdata function: ' obsdatafn])
  load(obsdatafn);
else
  disp(['loading obsdata'])
  [hydroobs,bathyobs,grid,waves8m,windEOP]=prepObsData(dnum,bathyfn,duck94Case);
  disp(['caching obsdata for next time: ' obsdatafn])
  save(obsdatafn,'hydroobs','bathyobs','grid','waves8m','windEOP');
end

% Phase-2 inversion code bathyAssim.m has been massaged into working with
% parfor, but is very memory-bound. (Hypothesis: parfor tries to allocate
% many copies of the 'bkgd' struct-array, which is about 4.5GB.)  Hence need
% to limit the number of workers, I seem to be able to get away with 2
% workers on plank which has 62GB RAM.
poolN=2;
currentPool=gcp('nocreate');
if(isempty(currentPool) | currentPool.NumWorkers ~= poolN)
  if(~isempty(currentPool))
    delete(gcp('nocreate'));
  end
  parpool('local',poolN);
end

% First phase-1 hydro-assimilating time loop, using default parameter values
modelinput=initModelInputs(duck94Case,grid,sedmodel);
bkgd=hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs);

% two-phase iterative assimlation
for iter=1:nitermax

  % phase 2, then phase 1
  [newparams,diagnostics]=bathyAssim(bkgd,bathyobs);
  modelinput.params=newparams;
  bkgd=hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs);

  % overview of results for this iteration
  clf, hold on
  cc='rbk';
  clear lstr
  for n=1:obsnt
    plot(grid.xFRF,bkgd2(bathyobs(n).obsn).hp,'-','color',cc(n))
    lstr{n}=['obs-time ' num2str(n)];
  end
  for n=1:obsnt
    plot(grid.xFRF,bkgd(bathyobs(n).obsn).hp,'--','color',cc(n))
  end
  for n=1:obsnt
    plot(grid.xFRF(bathyobs(n).h.ind),bathyobs(n).h.d,'o','color',cc(n))
  end
  title('OLD: dashed, NEW: solid, OBS: symbols')
  legend(lstr)
  set(gca,'ydir','r')

  % optional, save outputs
  if(dosave)
    outdir=['case_' duck94Case '_outputs/assimIter' num2str(iter)];
    unix(['mkdir -p ' outdir]);
    print('-dpng','-r300',[outdir '/bathyOutput.png'])
    bkgd_obsn = bkgd([bathyobs.obsn]);
    bkgd_1=bkgd(1);
    save([outdir '/output.mat'],'bkgd_obsn','bkgd_1','bathyobs','newparams','diagnostics')
  end

end  % two-phase iteration loop
