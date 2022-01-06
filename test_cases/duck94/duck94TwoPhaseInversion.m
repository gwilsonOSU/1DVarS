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
addpath(genpath('../../src'))  % hydroSedModel.m and its dependencies
addpath util  % helper functions
clearvars -except duck94Case sedmodel nitermax dosave parpoolN

%--------------------------------------
% user inputs
%--------------------------------------

% duck94 case: set this to a,b,c, or d
% duck94Case='d';
duck94Case

% sediment model: uncomment one of the following
if(~exist('sedmodel'))
  % sedmodel='dubarbier';
  sedmodel='vanderA';
end
sedmodel

% number of two-phase assimilation iterations
nitermax=8

% save figures and .mat files for each iteration
dosave=1;

% location of cache directories.  Each should have a ~5-10 GB available
% storage depending on run length
priorCacheDir='/tmp/bathyAssimCachePrior';
bkgdCacheDir='/tmp/bathyAssimCache';

% Parallel pool (used by bathyAssim.m). Note, if you want to go beyond ~12
% cores then the overhead of the parallel pool may bump up against available
% memory.  Some machines at OSU only have 25GB RAM and in that case you
% might run out of memory.
if(~exist('parpoolN'))
  parpoolN=12;
end
parpoolN
currentPool=gcp('nocreate');
if(isempty(currentPool) | currentPool.NumWorkers ~= parpoolN)
  if(~isempty(currentPool))
    delete(gcp('nocreate'));
  end
  parpool('local',parpoolN);
end

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

  % % OPTIONAL: case-c has a long quiescent period at the beginning.  For
  % % testing purposes, just skip to the storm period.
  % warning('saving time by starting case-c at oct3,0000EST. Be sure the obsdataCache reflects this!');
  % dnum(1)=datenum(1994,10,3,0,0,0);

 case 'd'

  warning('A rip channel developed on the bar after Oct 11, so case-d is suspect.')
  dnum(1)=datenum(1994,10,10,12,0,0);
  dnum(2)=datenum(1994,10,14,12,0,0);
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/1010_profile';

 case 'e'
  % Case-e is a concatenation of cases b-c, so covers the onshore and offshore
  % bar migration in sequence
  bathyfn='data/duck94_fulldataset/Dropbox/Duck94_profiles/0924_profile';
  dnum(1)=datenum(1994,9,24,12,0,0);
  dnum(2)=datenum(1994,10,4,12,0,0);

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

% redo bathyobs.h.e
warning('using obs-error 0.1m for bathy observations')
for n=1:length(bathyobs)
  bathyobs(n).h.e=ones(size(bathyobs(n).h.d))*.1;
end

% initial phase-1 hydro-assimilating time loop, using default parameter
% values.  Cache the results as the prior, which will be used throughout the
% iterations below.
modelinput=initModelInputs(duck94Case,grid,sedmodel);
hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs,priorCacheDir);

% two-phase iterative assimilation
for iter=1:nitermax
  disp(['starting iteration ' num2str(iter) ' of ' num2str(nitermax)])

  % For 1st iteration, prior==bkgd.  For subsequent iterations, prior~=bkgd.
  % Note, if you want to update the prior with each iteration, you can pass
  % thisPriorCache=='' for every iteration.  In that case, corrections are
  % allowed to stray far away from the initial prior, which may/not be
  % desirable (usually not).
  if(iter==1)
    thisBkgdCache=priorCacheDir;
  else
    thisBkgdCache=bkgdCacheDir;
  end

  % phase 2
  [newparams,diagnostics]=bathyAssim(bathyobs,priorCacheDir,thisBkgdCache);
  modelinput.params=newparams;
  modelinput.beta0=newparams.beta0;

  % if saving outputs, define the output directory
  if(dosave)
    outdir=['case_' duck94Case '_outputs/assimIter' num2str(iter)];
    unix(['mkdir -p ' outdir]);
  end

  % optional, save phase-2 results
  if(dosave)
    save([outdir '/phase2_output.mat'],'bathyobs','newparams','diagnostics')
    diary([outdir '/phase2_log.txt'])
    diagnostics.params0
    newparams
    diary off
  end

  % phase 1.  Note, this time we save the bkgd to bkgdCacheDir, not
  % priorCacheDir, so now bkgd~=prior
  hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs,bkgdCacheDir);

  % overview of results for this iteration
  clf, hold on
  cc='rbk';
  clear lstr
  obsnt=length(bathyobs);
  for n=1:obsnt
    bkgd0=load([priorCacheDir '/bkgd' num2str(bathyobs(n).obsn) '.mat']);
    plot(grid.xFRF,bkgd0.hp,'--','color',cc(n))
    lstr{n}=['obs-time ' num2str(n)];
  end
  for n=1:obsnt
    bkgd=load([bkgdCacheDir '/bkgd' num2str(bathyobs(n).obsn) '.mat']);
    plot(grid.xFRF,bkgd.hp,'-','color',cc(n))
  end
  for n=1:obsnt
    plot(grid.xFRF(bathyobs(n).h.ind),bathyobs(n).h.d,'o','color',cc(n))
  end
  title('PRIOR: dashed, NEW: solid, OBS: symbols')
  legend(lstr)
  set(gca,'ydir','r')
  xlim([100 400])

  % optional, save outputs
  if(dosave)
    print('-dpng','-r300',[outdir '/bathyOutput.png'])
    for n=1:length(bathyobs)
      bkgd_obsn(n) = load([bkgdCacheDir '/bkgd' num2str(bathyobs(n).obsn) '.mat']);
    end
    bkgd_1 = load([bkgdCacheDir '/bkgd1.mat']);
    save([outdir '/phase2_output.mat'],'bkgd_obsn','bkgd_1')
  end

end  % two-phase iteration loop

% clean up disk cache
disp('Deleting cached data')
unix(['rm ' priorCacheDir '/bkgd' num2str(n) '*.mat']);
unix(['rm ' bkgdCacheDir  '/bkgd' num2str(n) '*.mat']);
