%
% Single forward-run, used for manual testing
%
addpath util
addpath(genpath('../../src'))
clearvars -except duck94Case sedmodel

% USER-INPUT: Choose a case
% duck94Case='b';

if(~exist('sedmodel'))
  sedmodel='vanderA';
end

cacheDir='/tmp/bathyAssimCache';

%-------------------------------------------------
% End of user input.  No edits needed beyond this point.
%-------------------------------------------------

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

% load default model inputs
modelinput=initModelInputs(duck94Case,grid,sedmodel);

% make any desired modifications to model inputs
modelinput.params.Cf=0;  % turn off slope-based suspended transport

% run forward model (with hydro assimilation)
numsubsteps=1;
doplot=1;
hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs,cacheDir,numsubsteps,doplot);

% overlay final bathy
subplot(321)
hold on
plot(grid.xFRF(bathyobs(end).measind),bathyobs(end).h.d,'k-','linewidth',1.5)
print -dpng casee_defaults_WBLSplittingMethod_noassim.png

disp('continue manually')
return;

% load cached data into memory and delete from disk
disp('loading cache files')
nt=length(hydroobs)-1;
for n=1:nt
  fn=[cacheDir '/bkgd' num2str(n) '.mat'];
  bkgd(n)=load(fn);
  unix(['rm -f ' fn]);
end
