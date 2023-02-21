% Prediction of transport for Silva et al. (2011), cases with waves and
% current combined
%
addpath vanderA2013
clear

%-------------------------------------
% tuning params for van der A et al. (2013)
%-------------------------------------

params.n=1.2;  % default 1.2.  Larger values promote offshore bar migration
params.m=11;  % default 11; hsu et al. 11.  Just scales everything up
params.xi=0;  % for progressive waves only, O(1) according to Kranenburg (2013).  VDA pg. 35 says xi=1.7
params.alpha=8.2;  % default 8.2.  Comes in eqn 27-28, not the same as eqn 19

%-------------------------------------
% experimental params from Silva et al. (2011)
%-------------------------------------

% constants
d50 =0.2/1000;  % median grainsize, m
sigg=1.18;  % geometric stdev for lognormal sediment distribution

% params from Table 1 of paper
U0 = [0.00 0.00 0.00 0.00 -0.22 -0.44 -0.22 -0.44 0.00 0.00 0.00];  % "undertow", m/s, -'ve valued
urms=[90 88 88 86 89 88 86 86 86 94 87]/100;  % rms wave velocity, m/s
T=[7 10 7 10 7 7 7 7 7 10 7];  % wave period, s
beta=[64 63 72 72 64 64 71 71 61 60 53]/100;  % accel skewness, nondim
R=[50 49 50 49 50 50 50 51 59 59 59];  % vel skewness nondim
N=[5 2 5 4 7 3 4 3 4 2 4];  % nondim

% measured transport rates (Table 1)
qs_lab = [ +.0539 +.0443 +.1137 +.0847 -.1133 -.3826 -.0489 -.2241 +.1845 +.2797 +.1244];  % kg/m/s
rhop = 2650;  % kg/m3
qs_lab = qs_lab / rhop;  % m2/s

%-------------------------------------
% derived params
%-------------------------------------

ncases=length(U0);

% mean current, 2d version
udelta(:,1)=U0;
udelta(:,2)=zeros(ncases,1);

% d90, assuming lognormal distribution
sig=log(sigg);
mu=log(d50);
d90 = exp( mu + sqrt(2)*sig*erf(2*.9-1) );

% settling velocity
addpath ../../hydroSedModel/ws_brownLawler2003
ws=ws_brownLawler(.8*d50);
wsc=ws;  % constant ws, no effect of orbital vels
wst=ws;  % constant ws, no effect of orbital vels

% wave shape parameters as defined by VDA2013.  Code copied from my
% parameterized version of qtrans_vanderA.m
r_r2012 = 2*beta-1;  % eqn 5
phi_r2012 = -pi/4;  % see end of section 2.2
omega=2*pi./T;
phiuc=asin(-r_r2012.*sin(phi_r2012)./(1+sqrt(1-r_r2012.^2)))./omega;   % phase at first upcrossing
phidc=pi-phiuc;  % phase at first downcrossing
tuc=phiuc./omega;  % time of first upcrossing
tdc=phidc./omega;  % time of first downcrossing
for i=1:ncases
  tcr(i)=Uwave_ruessink2012_tcrest(omega(i),r_r2012(i),phi_r2012,T(i)*1/4);  % time of crest
  ttr(i)=Uwave_ruessink2012_tcrest(omega(i),r_r2012(i),phi_r2012,T(i)*3/4);  % time of trough
end
Tc = tdc-tuc;  % duration of crest
Tt=T-Tc;     % duration of trough
Ttu=ttr-tdc;  % duration of deceleration under trough
Tcu=tcr-tuc;  % duration of acceleration under crest

% Abreu et al. (2010) wave velocity time series formula
doplot=0;  % optional, visually verify the crest/trough timings gotten above
phs=linspace(0,2*pi,1000);
for i=1:ncases

  % calculate u(t) from Abreu equations
  f=sqrt(1-r_r2012(i).^2);
  f1 = sin(phs) + r_r2012(i).*sin(phi_r2012)./(1+sqrt(1-r_r2012(i).^2));
  f2 = 1 - r_r2012(i).*cos(phs+phi_r2012);
  u = f.*f1./f2;
  u = u*urms(i)/sqrt(mean(u.^2));

  % get crest and trough velocities
  uhat(i)=sqrt(2*mean(u.^2));   % rms wave velocity for full wave cycle, eqn 8
  uhatc(i)=+max(u);  % max velocity under crest, see text below eqn 7
  uhatt(i)=-min(u);  % max velocity under trough (positive valued)

  % optional, show the results
  if(doplot)
    plot(phs/omega(i),u)
    hold on
    ax=axis;
    plot(tuc(i)*[1 1],ax(3:4),'k--')
    plot((tuc(i)+Tc(i))*[1 1],ax(3:4),'k--')
    plot((tuc(i)+Tc(i)+Ttu(i))*[1 1],ax(3:4),'k--')
    plot((tuc(i)+Tcu(i))*[1 1],ax(3:4),'k--')
    hold off
    drawAxis;
    pause;
  end

end

%-------------------------------------
% forward model analysis
%-------------------------------------

for i=1:ncases
  [qs(i),bkgd(i)]=qtrans_vanderA(d50,d90,wsc,wst,udelta(i,:),uhat(i),uhatc(i),uhatt(i),T(i),Tc(i),Tt(i),Ttu(i),Tcu(i),999,params);
end

% reproduce VDA fig. 8.  Show the cases with zero current (U0==0) as gray
clf, hold on
ind=find(abs(U0)==0)
plot(qs_lab(ind)*1e6,qs(ind)*1e6,'p','color',[1 1 1]*.6)
ind=find(abs(U0)~=0)
plot(qs_lab(ind)*1e6,qs(ind)*1e6,'pk')
axis equal
axis([-1 1 -1 1]*250)
drawAxis
hold on
ax=axis;
plot([-1 1]*250,   [-1 1]*250,'k')
plot([-1 1]*250,2.*[-1 1]*250,'k--')
plot([-1 1]*250,.5*[-1 1]*250,'k--')
axis(ax)
box on
set(gca,'xtick',-250:50:250)
set(gca,'ytick',-250:50:250)
xlabel('Measured (x10^6 m^2/s)')
ylabel('Predicted (x10^6 m^2/s)')
legend('Silva (2011), no current',...
       'Silva (2011), opposing current',...
       'location','northwest')
