function modelinput=initModelInputs(duck94Case,grid,sedmodel)
%
% modelinput=initModelInputs(duck94Case,grid,sedmodel)
%
% Initializes modelinput struct required for hydro-assimilation time-loop.
% Parameter values will be set to defaults.  If you like you can over-ride
% these default values by overwriting them before proceeding to the time
% loop.

% initialize sediment transport params
params=struct;
params.fv=0.101;  % breaking-induced eddy viscosity calib parameter, default 0.101
params.ks=0.0082; %0.0082;  % roughness calib parameter, default 0.0082 m
params.lambda=1.57;  % Dubarbier et al.'s eqn 9 to shift velocities shoreward, default 1.57
if(strcmp(sedmodel,'dubarbier'))
  params.Cw=.002;  % default 0.000483;  % hsu et al. 0.0046
  params.Cc=0.01;  % default 0.02002;  % hsu et al. 0.0053
  params.Cf=0.01;  % default 0.01173 ; hsu et al. 0.0
  params.Ka=.75e-4;  % default 0.631e-4;  % hoefel & elgar 1.4e-4; hsu et al. 0.0
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;  % vanderA 1.2.  Larger values -> offshore transport
  params.m=11; %11;  % vanderA 11; hsu et al. 11.  Just scales everything up
  params.xi=1.7;  % O(1) according to Kranenburg (2013).  VDA pg. 35 says 1.7
  params.alpha=8.2;  % comes in eqn 27-28, not the same as eqn 19.  Default 8.2

  % NOTE: if Cc and Cf are omitted, then qtrans_vanderA.m will calculate above-WBL transport automatically (recommended)
  % params.Cc=0.005;  % suspended sediment stirring+undertow effect.  Default 0.01
  % params.Cf=0.01;  % suspended sediment stirring+slope effect. Default 0.01, but 0.03 may be good

  params.streamingType='v';  % 'n' for nielsen2006 BL-streaming, 'v' for VDA13's formulation, or '0' for off
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% % OLD: If not using suspended sediment in vanderA model, will need to
% % amplify the undertow to get sufficient offshore transport for the offshore
% % bar migration events.  The below parameters do the trick.
% if(duck94Case~='b')
%   warning('using tuned undertow for offshore migration')
%   params.fv=.03;
%   params.ks=2.5*180e-6;  % 2.5*d50
% end

% initialize other model params
modelinput=grid;  % initialize: x,h,xFRF
modelinput.ka_drag=0.0125;  % tuned by Ruessink et al. (2001)
modelinput.beta0=0.05;
modelinput.d50=zeros(grid.nx,1);
modelinput.d90=zeros(grid.nx,1);
modelinput.d50(grid.xFRF<150)=400e-6;  % shoreface coarse sand, Birkemeier et al. (1985)
modelinput.d90(grid.xFRF<150)=400e-6;  % shoreface coarse sand
modelinput.d50(grid.xFRF>=150)=180e-6; % Dubarbier et al. (2015) used 170um here, but 180 is from Birkemeier et al. (1985)
modelinput.d90(grid.xFRF>=150)=240e-6;  % Birkemeier et al., 1985
modelinput.params=params;
modelinput.sedmodel=sedmodel;
modelinput.betaType='const';
modelinput.gammaType=2003;

% initialize covariances for hydro-assimilation step
xx=meshgrid(grid.xFRF);
sig_h=0*exp(-3*(grid.xFRF-200).^2/150^2);  % set to zero, no h-corrections in hydro assimlation step
Ch0=diag(sig_h)*exp(-3*(xx-xx').^2/50^2)*diag(sig_h);
modelinput.Ch=Ch0*0;  % set to zero, no h-corrections in hydro assimlation step
modelinput.Cs=0;  % set to zero, don't correct sed params in hydro assim step
modelinput.Cgamma=.1^2*exp(-3*(xx-xx').^2/25^2);  % correct wave breaking errors
modelinput.CdAw=.1^2*exp(-3*(xx-xx').^2/50^2);  % correct wave asymmetry errors
modelinput.CdSw=.1^2*exp(-3*(xx-xx').^2/50^2);  % correct wave skewness errors
modelinput.CH0=0.25^2;  % correct wave BC errors
modelinput.Ctheta0=deg2rad(10)^2;  % correct wave BC errors
modelinput.Cka=0.005^2;  % correct bottom friction (only relevant if assimilating vbar)

% % TEST: disable hydro assimilation
% warning('setting all hydro covariances to zero!!');
% modelinput.Ch=Ch0*0;  % set to zero, no h-corrections in hydro assimlation step
% modelinput.Cs=0;  % set to zero, don't correct sed params in hydro assim step
% modelinput.Cgamma=0*.1^2*exp(-3*(xx-xx').^2/25^2);  % correct wave breaking errors
% modelinput.CdAw=0*.1^2*exp(-3*(xx-xx').^2/50^2);  % correct wave asymmetry errors
% modelinput.CdSw=0*.1^2*exp(-3*(xx-xx').^2/50^2);  % correct wave skewness errors
% modelinput.CH0=0*0.25^2;  % correct wave BC errors
% modelinput.Ctheta0=0*deg2rad(10)^2;  % correct wave BC errors
% modelinput.Cka=0*0.005^2;  % correct bottom friction (only relevant if assimilating vbar)
