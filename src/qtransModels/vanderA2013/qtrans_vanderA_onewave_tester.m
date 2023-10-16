clear

addpath ../..  % need to add the path that contains physicalConstants.m

% Values for wave parameters.  These are somewhat arbitrarily assigned here
% just to get something in the right ballpark for the purposes of testing
% the code.  Normally you would get these from measurements.
uw=sin(linspace(0,2*pi,1000));
uhat=sqrt(2*mean(uw.^2));
uhatc=max(uw);
uhatt=-min(uw);
T=10;
Tc=5;
Tt=5;
Tcu=2;
Ttu=2;
c=4;
d50=180e-6;
d90=400e-6;
h=1;
tanbeta=0.001;
udelta=[0.1 0];
wsc=0.01;
wst=0.01;

% inclusion of the "previous wave" Omegatc is optional.  Would normally get
% this from previous wave's calculation.  If Omegatc is not specified, the
% model will assume the waves are periodic, i.e. the current wave's Omegatc
% is used to represent the previous wave (as in VDA13).
Omegatc=1e-4;

% params struct, here are the defaults...
param.n=1.2;  % vanderA 1.2.  Larger values -> offshore transport
param.m=11;   % vanderA 11; hsu et al. 11.  Just scales everything up
param.xi=1.7;  % O(1) according to Kranenburg (2013).  VDA pg. 35 says 1.7
param.alpha=8.2;  % comes in eqn 27-28, not the same as eqn 19.  Default 8.2
param.streamingType='v';  % 'n' for nielsen2006 BL-streaming, 'v' for VDA13's formulation (default), or '0' for off

% If param.Cc and param.Cf are omitted, then the code will calculate
% above-WBL transport automatically, with default parameter values
% (recommended).  Comment them out below to get this default behavior.  Or,
% you can uncomment them and tweak the values to try and improve the
% results.
%
% params.Cc=0.005;  % suspended sediment stirring+undertow effect.  Default 0.01
% params.Cf=0.01;  % suspended sediment stirring+slope effect. Default 0.01

% Yet another option is to fully disable the calculation of above-WBL
% transport.  Set this by defining param.nosusp=1.  The default is
% param.nosusp=0 (i.e., include above-WBL transport), and this default is
% chosen if you don't define the nosusp field in your param struct.
param.nosusp=0;

% The undertow-driven transport results may be sensitive to the value of the
% boundary layer thickness delta, since this determines how much sediment is
% assigned to be transported inside vs. above the WBL.  VDA's paper just
% sets it as a constant, delta=0.2m.  A possibility would be to run
% udelta_reniers2004.m to get a more-justified (model-based) estimate of
% delta, or try to estimate it from measurements.
delta=0.2;

% run the onewave code (note, optionally uncomment the input Omegatc to
% include it in the calculation).
[qs,workspc,OmegatcOut]=qtrans_vanderA_onewave(uw,uhat,uhatc,uhatt,T,Tc,Tt,Tcu,Ttu,c,d50,d90,h,tanbeta,udelta,delta,wsc,wst,param); %,Omegatc);
