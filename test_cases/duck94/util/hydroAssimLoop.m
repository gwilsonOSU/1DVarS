function bkgdall=hydroAssimLoop(modelinput,grid,waves8m,windEOP,obs)
%
% bkgd=hydroAssimLoop(modelinput,grid,waves8m,windEOP,hydroobs)
%
% Phase-1 assimilation.  Run a time-loop to predict bathymetry, while
% assimilating hydro data to keep wave model errors in check.
%
% INPUT:
%
% modelinput: struct containing initial model inputs, e.g. see initModelInputs.m
% waves8m: 8m-array wave data, see prepObsData.m
% windEOP: wind data from end of FRF pier, see prepObsData.m
% hydroobs: time-dependent observation data to be assimilated, see prepObsData.m
%
% OUTPUT:
%
% bkgd: An array of structs with model outputs, one struct for each
% observation time step.  Its fields are described in more detail in
% hydroAssimOneStep.m and hydroSedModel.m.
%

% hard-coded options
doCovUpdate=0;  % not updating Ch over time, so set this to 0
nsubsteps=1;  % number of time-sub-steps for hydroSedModel.m.  Sometimes
              % improves stability. but =1 seems to work fine for duck94
              % tests

nx=grid.nx;

% begin time loop
for n=1:(length(obs)-1)
  disp(['time ' num2str(n) ' of ' num2str(length(obs))])

  % Set wave conditions for this time step.  Note, use Tm1 wave period
  % (corresponding to waves8m.sigmam) following Ruessink et al. (2012)
  % formulation for wave shape
  modelinput.H0    =interp1(waves8m.dnum_est,waves8m.Hrms  ,obs(n).dnum_est);
  modelinput.theta0=interp1(waves8m.dnum_est,waves8m.theta0,obs(n).dnum_est);
  modelinput.omega =interp1(waves8m.dnum_est,waves8m.sigmam,obs(n).dnum_est);
  modelinput.tau_wind=interp1(windEOP.dnum_est,windEOP.tau,obs(n).dnum_est);
  modelinput.tau_wind=repmat(modelinput.tau_wind,nx,1);
  modelinput.dgamma=zeros(nx,1);
  modelinput.dAw=zeros(nx,1);
  modelinput.dSw=zeros(nx,1);

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
    detady = 1/mean(g*thish)*( Fy - modelinput.tau_wind(1,2)/rho );  % solve longshore momentum budget
    modelinput.detady=detady*ones(nx,1);
  end

  % Run model to step forward in time, and assimilate data to correct hydro
  % errors.  This takes the initial 'modelinput' struct, and produces a new
  % 'bkgd' struct containing a forecast valid for the current time step,
  % along with updated covariances for this time step.
  dt=diff([obs(n+[0:1]).dnum_est])*24*60*60;
  % verb=1; figure(1);
  verb=0;
  bkgd=hydroAssimOneStep(modelinput,thisobs,dt,nsubsteps,verb);
  bkgd.dnum_est=obs(n).dnum_est;
  bkgd.tide=tide;
  bkgd.h =bkgd.h -tide;
  bkgd.hp=bkgd.hp-tide;
  bkgdall(n)=bkgd;

  % Use the forecasted bathymetry as inputs for the next time step t(n+1).
  % Note, this code does not attempt to propagate bathymetry covariance
  % since it is not correcting bathymetry.  Therefore Ch stays the same.
  modelinput.h=bkgd.hp(:,end);

  %-----------------------------------------
  % plot results for this time step
  %-----------------------------------------

  % figure(2)
  clf
  lw=1.5;
  subplot(321), hold on
  plot(grid.xFRF,grid.h,'g','linewidth',lw)
  plot(grid.xFRF,bkgd.h,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).h.ind),obs(n).h.d,'ko')
  set(gca,'ydir','r')
  ylabel('h [m]')
  subplot(322), hold on
  plot(grid.xFRF,bkgd.Hrms,'b','linewidth',lw)
  plot(grid.xFRF(obs(n).H.ind),obs(n).H.d,'ko')
  ylabel('H_{rms} [m]')
  subplot(323), hold on
  plot(grid.xFRF,bkgd.Q,'b','linewidth',lw)
  ylabel('Q [m^2/s]')
  ylim([-1 1]*1e-4)
  subplot(324), hold on
  plot(grid.xFRF,bkgd.ubar(:,2),'b','linewidth',lw)
  plot(grid.xFRF,bkgd.udelta(:,2),'b--','linewidth',lw)
  plot(grid.xFRF(obs(n).v.ind),obs(n).v.d,'ko')
  ylabel('v [m/s]')
  subplot(325), hold on
  plot(grid.xFRF,bkgd.ubar(:,1),'b','linewidth',lw)
  plot(grid.xFRF,bkgd.udelta(:,1),'b--','linewidth',lw)
  plot(grid.xFRF(obs(n).u.ind),obs(n).u.d,'ko')
  ylim([-.5 .1])
  ylabel('u [m/s]')
  subplot(326), hold on
  plot(grid.xFRF,bkgd.Aw(:,1),'b','linewidth',lw)
  plot(grid.xFRF(obs(n).A.ind),obs(n).A.d,'bo')
  plot(grid.xFRF,bkgd.Sw(:,1),'r','linewidth',lw)
  plot(grid.xFRF(obs(n).S.ind),obs(n).S.d,'ro')
  ylim([-1 1])
  ylabel('A_w,S_w')
  for j=1:6
    subplot(3,2,j)
    box on
    xlim([100 400])
    ax=axis;
    if(j==1)
      [~,ishore]=min(abs(bkgd.h));
      plot(ax(1:2),-[1 1]*0,'k--')
    else
      plot(ax(1:2),[0 0],'k--')
    end
    plot(grid.xFRF(ishore)*[1 1],ax(3:4),'k--')
    axis(ax)
  end
  subplot(321)
  title(['yday ' num2str(obs(n).dnum_est-datenum(1994,1,0)) ...
        ', ' datestr(obs(n).dnum_est) 'EST'])
  legend('initial, n=1',['n=' num2str(n) ' of ' num2str(length(obs))],'observed')
  pause(.01)

end
