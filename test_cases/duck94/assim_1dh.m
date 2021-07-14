function [posterior,fcst]=assim_1dh(prior,obs,dt,verb,doCovUpdate)

if(~exist('verb'))
  verb=1;
end
if(~exist('doCovUpdate'))
  doCovUpdate=1;
end

nitermax=500;

% unpack some variables for convenience
Chs=prior.Chs;
Cgamma=prior.Cgamma;
CH0=prior.CH0;
Ctheta0=prior.Ctheta0;
Cka=prior.Cka;
x=prior.x;
nx=length(prior.x);
ns=size(Chs,1)-nx;  % number of sediment params
Ch=Chs(1:nx,1:nx);

% if any obs types missing, set to empty
fld='Hvh';
noobs.ind=[];
noobs.d=[];
noobs.e=[];
for i=1:length(fld)
  if(~isfield(obs,fld(i)))
    obs=setfield(obs,fld(i),noobs);
  end
end

% initialize bkgd hydro model state using the prior.  Note, assume dgamma=0 in prior.
bkgd=hydro_ruessink2001(prior.x,prior.h,...
                        prior.H0,prior.theta0,prior.omega,prior.ka_drag,...
                        prior.tauw,prior.detady);
horig=prior.h;
for fld=fields(bkgd)'
  fld=cell2mat(fld);
  prior=setfield(prior,fld,getfield(bkgd,fld));
end
prior.h=horig;  % don't cut off dry points

% outer loop
eps=nan;
for n=1:nitermax
  disp(['iteration ' num2str(n) ', itermax = ' num2str(nitermax) ', eps = ' num2str(eps)])

  % run forecast for this outer loop iteration, linearized about most-recent
  % background state.  The latter (bkgd) is updated witht each outer loop
  % iteration.
  fcst=prior;
  tl_h=fcst.h-bkgd.h;
  tl_H0=fcst.H0-bkgd.H0;
  tl_theta0=fcst.theta0-bkgd.theta0;
  tl_ka_drag=fcst.ka_drag-bkgd.ka_drag;
  tl_omega=0; % no uncertainty in omega
  tl_tauw=0*bkgd.tauw;  % no uncertainty in tauw
  tl_detady=0*bkgd.detady; % no uncertainty in detady
  tl_dgamma=-bkgd.dgamma; % note, assume prior dgamma==0
  [tl_H,tl_theta,tl_v,tl_k,tl_Ew,tl_Er,tl_Dr] = ...
      tl_hydro_ruessink2001(tl_h,tl_H0,tl_theta0,tl_omega,...
                            tl_ka_drag,tl_tauw,tl_detady,tl_dgamma,bkgd);
  fcst.Hrms=bkgd.Hrms+tl_H;
  fcst.vbar=bkgd.vbar+tl_v;
  fcst.theta=bkgd.theta+tl_theta;
  fcst.k=bkgd.k+tl_k;
  fcst.Ew=bkgd.Ew+tl_Ew;
  fcst.Er=bkgd.Er+tl_Er;
  fcst.Dr=bkgd.Dr+tl_Dr;

  % initialize representer sub-matrices
  %
  % LEGEND:
  %
  % R_XY = sensitivity of observation type Y, to delta-perturbations of
  % observation X.  These are sub-blocks of L*M*C*M'*L' (note, C is the
  % prior covariance).
  %
  % r_X = sensitivity of model input vector phi = [h; sedparams; H0; theta0; ka_drag])
  % to delta-perturbations of observation X.  These are rows of C*M'*L'.
  %
  clear R_* r_*
  r_H=zeros(2*nx+ns+3,length(obs.H.ind));
  r_v=zeros(2*nx+ns+3,length(obs.v.ind));
  r_h=zeros(2*nx+ns+3,length(obs.h.ind));
  R_HH=zeros(length(obs.H.ind),length(obs.H.ind));
  R_Hv=zeros(length(obs.H.ind),length(obs.v.ind));
  R_Hh=zeros(length(obs.H.ind),length(obs.h.ind));
  R_vH=zeros(length(obs.v.ind),length(obs.H.ind));
  R_vv=zeros(length(obs.v.ind),length(obs.v.ind));
  R_vh=zeros(length(obs.v.ind),length(obs.h.ind));
  R_hH=zeros(length(obs.h.ind),length(obs.H.ind));
  R_hv=zeros(length(obs.h.ind),length(obs.v.ind));
  R_hh=zeros(length(obs.h.ind),length(obs.h.ind));

  % compute representers for wave height: apply delta-perturbations to
  % observations of type X, to compute (a) sensitivity of model inputs
  % (r_* matrices), and (b) sensitivity of observations of type Y (R_*
  % matrices)
  for i=1:length(obs.H.ind)
    comb=zeros(nx,1);
    comb(obs.H.ind(i))=1;  % data functional (delta-fn, aka identity matrix)

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauw,ad_detady,ad_dgamma] = ...
        ad_hydro_ruessink2001(comb,0*comb,0*comb,0*comb,0*comb,0*comb,0*comb,bkgd);

    % multiply into covariances to get the matrix of representers C*M'
    r_H(:,i)=[Chs*[ad_h; zeros(ns,1)]; Cgamma*ad_dgamma; CH0*ad_H0; Ctheta0*ad_theta0; Cka*ad_ka_drag];

    % apply TL model and measure it to get representer matrix, M*C*M'
    [tl_H,tl_theta,tl_v]=tl_hydro_ruessink2001(Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,0*ad_omega,...
                                               Cka*ad_ka_drag,0*ad_tauw,0*ad_detady,Cgamma*ad_dgamma,bkgd);
    R_HH(i,:)=tl_H(obs.H.ind);
    R_Hv(i,:)=tl_v(obs.v.ind);
    R_Hh(i,:)=r_H(obs.h.ind,i);

  end

  % compute representers for longshore current, same operations as for wave
  % height in previous block of code
  for i=1:length(obs.v.ind)
    comb=zeros(nx,1);
    comb(obs.v.ind(i))=1;  % data functional (delta-fn, aka identity matrix)

    % I*M', where M is the TL model and I is identity
    [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_tauw,ad_detady,ad_dgamma] = ...
        ad_hydro_ruessink2001(0*comb,0*comb,comb,0*comb,0*comb,0*comb,0*comb,bkgd);

    % multiply into covariances to get the matrix of representers C*M'
    r_v(:,i)=[Chs*[ad_h; zeros(ns,1)]; Cgamma*ad_dgamma; CH0*ad_H0; Ctheta0*ad_theta0; Cka*ad_ka_drag];

    % apply TL model and measure it to get representer matrix, M*C*M'
    [tl_H,tl_theta,tl_v]=tl_hydro_ruessink2001(Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,0*ad_omega,...
                                               Cka*ad_ka_drag,0*ad_tauw,0*ad_detady,Cgamma*ad_dgamma,bkgd);
    R_vH(i,:)=tl_H(obs.H.ind);
    R_vv(i,:)=tl_v(obs.v.ind);
    R_vh(i,:)=r_v(obs.h.ind,i);

  end

  % compute representers for h, same operations as in previous block(s) of code
  for i=1:length(obs.h.ind)
    comb=zeros(nx,1);
    comb(obs.h.ind(i))=1;  % data functional (delta-fn, aka identity matrix)

    % I*M', where M is the TL model and I is identity.  No need to run a model
    % since h is a static parameter in the model
    ad_h=comb;
    ad_H0=0;
    ad_theta0=0;
    ad_omega=0;
    ad_ka_drag=0;
    ad_tauw=zeros(nx,2);
    ad_detady=zeros(nx,1);
    ad_dgamma=zeros(nx,1);

    % multiply into covariances to get the matrix of representers C*M'
    r_h(:,i)=[Chs*[ad_h; zeros(ns,1)]; Cgamma*ad_dgamma; CH0*ad_H0; Ctheta0*ad_theta0; Cka*ad_ka_drag];

    % apply TL model and measure it to get representer matrix, M*C*M'
    [tl_H,tl_theta,tl_v]=tl_hydro_ruessink2001(Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,0*ad_omega,...
                                               Cka*ad_ka_drag,0*ad_tauw,0*ad_detady,Cgamma*ad_dgamma,bkgd);
    R_hH(i,:)=tl_H(obs.H.ind);
    R_hv(i,:)=tl_v(obs.v.ind);
    R_hh(i,:)=r_h(obs.h.ind,i);

  end

  % assemble matrices for updating
  Cd=diag([obs.H.e(:);
           obs.v.e(:);
           obs.h.e(:)].^2);
  CMt=[r_H r_v r_h];
  d=[obs.H.d(:);
     obs.v.d(:);
     obs.h.d(:)];
  Lu=[fcst.Hrms(obs.H.ind);
      fcst.vbar(obs.v.ind);
      fcst.h(obs.h.ind)];
  R=[R_HH R_Hv R_Hh;
     R_vH R_vv R_vh;
     R_hH R_hv R_hh];
  obstype=[repmat('H',[length(obs.H.d) 1]);
           repmat('v',[length(obs.v.d) 1])
           repmat('h',[length(obs.h.d) 1])];

  % compute the update.  Ignore sediment params in this update
  hedge=1; %tanh(n/5);  % reduce magnitude of update for 1st few iterations
  update=hedge*CMt*inv(R+Cd)*(d-Lu);
  if(n>1)
    update=.5*(u0+update);
  end
  u0=update;
  posterior=prior;
  posterior.h = prior.h + update(0*nx+[1:nx]);
  % note, sed param updates would be update(1*nx+[1:ns])
  posterior.dgamma = prior.dgamma + update(1*nx+ns+[1:nx]);
  posterior.H0=prior.H0+update(2*nx+ns+1);
  posterior.theta0=prior.theta0+update(2*nx+ns+2);
  posterior.ka_drag=prior.ka_drag+update(2*nx+ns+3);
  posterior.ka_drag=max(1e-4,posterior.ka_drag);  % limiter on ka_drag
  bkgd0=bkgd;

  % obtain new background state
  bkgd=hydro_ruessink2001(posterior.x,posterior.h,...
                          posterior.H0,posterior.theta0,...
                          posterior.omega,posterior.ka_drag,...
                          posterior.tauw,posterior.detady,posterior.dgamma);

  % show the results
  if(verb)
    clf
    subplot(221), hold on
    plot(prior.x,prior.h,'r')
    plot(prior.x,prior.h-sqrt(diag(prior.Chs(1:nx,1:nx))),'r--')
    plot(prior.x,prior.h+sqrt(diag(prior.Chs(1:nx,1:nx))),'r--')
    plot(posterior.x,posterior.h,'b')
    % plot(posterior.x,posterior.h-posterior.hErr,'b--')
    % plot(posterior.x,posterior.h+posterior.hErr,'b--')
    set(gca,'ydir','reverse')
    ylabel('h [m]')
    % legend('prior','inverted')
    subplot(222), hold on
    plot(prior.x,prior.Hrms,'r')
    plot(posterior.x,bkgd.Hrms,'b')
    errorbar(x(obs.H.ind),obs.H.d,obs.H.e,'ko')
    ylabel('H_{rms} [m]')
    title(['H0 prior=' num2str(prior.H0,2) 'm, final=' num2str(posterior.H0,2) 'm']);
    subplot(223), hold on
    plot(prior.x,rad2deg(prior.theta),'r')
    plot(posterior.x,rad2deg(bkgd.theta),'b')
    ylabel('theta [deg]')
    title(['theta0 prior=' num2str(rad2deg(prior.theta0),2) 'deg, final=' num2str(rad2deg(posterior.theta0),2) 'deg']);
    subplot(224), hold on
    plot(prior.x,prior.vbar,'r')
    plot(posterior.x,bkgd.vbar,'b')
    errorbar(x(obs.v.ind),obs.v.d,obs.v.e,'ko')
    ylabel('v [m/s]')
    title(['ka prior=' num2str(prior.ka_drag,3) 'm, final=' num2str(posterior.ka_drag,3) 'm']);
    for i=1:4
      subplot(2,2,i)
      box on
      % drawAxis
    end
    pause(.1)
  end

  % check for convergence
  eps=sum((bkgd0.h-bkgd.h).^2)/trace(prior.Chs(1:nx,1:nx));
  if(CH0>0)
    eps=eps+(bkgd0.H0-bkgd.H0)^2/prior.CH0;
  end
  if(Ctheta0>0)
    eps=eps+(bkgd0.theta0-bkgd.theta0)^2/prior.Ctheta0;
  end
  if(Cka>0)
    eps=eps+(bkgd0.ka_drag-bkgd.ka_drag)^2/prior.Cka;
  end
  if(eps<1e-5)  % TODO decrease to 1e-4 for final version
    break;
  end

end  % outer loop iterations

% update sediment parameters from the final outer-loop NL model run
posterior0=posterior;
for fld=fields(bkgd)'
  fld=cell2mat(fld);
  posterior=setfield(posterior,fld,getfield(bkgd,fld));
end
posterior.h=posterior0.h;  % bkgd version has hmin cutoff
fld=fields(posterior.params);
for i=1:length(fld)
  posterior.params=setfield(posterior.params,fld{i},getfield(posterior.params,fld{i})+update(nx+i));
end

% calculate posterior covariances
C2=blkdiag(Chs,Cgamma,CH0,Ctheta0,Cka)-CMt*inv(R+Cd)*CMt';
posterior.Chs=C2(1:(nx+ns),1:(nx+ns));
posterior.Cgamma=C2(1*nx+ns+[1:nx],1*nx+ns+[1:nx]);
posterior.CH0    =C2(2*nx+ns+1);
posterior.Ctheta0=C2(2*nx+ns+2);
posterior.Cka    =C2(2*nx+ns+3);

% forecast h for the next obs-time t+dt
disp('running deterministic forecast')
fcst = hydroSedModel(posterior.x,posterior.h,...
                     posterior.H0,posterior.theta0,posterior.omega,...
                     posterior.ka_drag,posterior.tauw,posterior.detady,...
                     posterior.dgamma,...
                     posterior.d50,posterior.d90,posterior.params,posterior.sedmodel,dt);
posterior.hp=fcst.hp;
posterior.Q =fcst.Q ;

% calculate covariance of the forecast bathymetry, and its covariance with
% sediment transport parameters.  This will facilitate parameter corrections
% in future analysis cycles
if(~doCovUpdate)  % skip this step
  disp('not calculating forecast covariance Chsp')
  return;
end
disp('propagating forecast covariance')
Chsp=nan(nx+ns);

parfor i=1:nx
  % disp(['  gridpoint ' num2str(i) ' of ' num2str(nx)])
  if(floor(i/nx*10)>floor((i-1)/nx*10))
    disp([num2str(floor(i/nx*10)*10) '%'])
  end

  % adjoint model acting on identity matrix
  comb=zeros(nx,1);
  comb(i)=1;  % data functional (delta-fn, aka identity matrix)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_dgamma, ...
      ad_tau_wind,ad_detady,ad_d50,ad_d90,ad_params] = ...
      ad_hydroSedModel(0*comb,0*comb,0*comb,0*comb,0*comb,comb,fcst);
 
  % convert sediment transport params from struct to vector
  fld=fields(ad_params);
  ad_pp=zeros(length(fld),1);     % Old line

  for j=1:length(fld)
    ad_pp(j)=getfield(ad_params,fld{j});  
  end

  % apply the full posterior covariance matrix to the adjoint vector, then
  % re-package the outputs for input to TL model
  phi=C2*[ad_h; ad_pp; ad_dgamma; ad_H0; ad_theta0; ad_ka_drag];
  ad_h=phi(1:nx);
  for j=1:length(fld)
    ad_params=setfield(ad_params,fld{j},phi(nx+j));
  end
  ad_dgamma = phi((nx+ns+1):(2*nx+ns));   
  ad_H0=phi(nx+ns+nx+1);
  ad_theta0=phi(nx+ns+nx+2);
  ad_ka_drag=phi(nx+ns+nx+3);

  % tangent linear run, I*C*I'
  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp] = ...
      tl_hydroSedModel(ad_h,...
                       ad_H0,...
                       ad_theta0,...
                       0*ad_omega,...
                       ad_ka_drag,...
                       ad_dgamma,...
                       0*ad_tau_wind,...
                       0*ad_detady,...
                       0*ad_d50,...
                       0*ad_d90,...
                       ad_params,...
                       fcst);

  % extract the forecast bathy covariance for this gridpoint
  Chsp(:,i)=[tl_hp; phi(nx+[1:ns])];

end  % loop over gridpoints
parfor i=1:ns

  % note, an adjoint model run is not required here since the sediment
  % transport params are an input not an output
  ad_h=zeros(nx,1);
  ad_dgamma = zeros(nx,1);
  ad_H0=0;
  ad_theta0=0;
  ad_omega=0;
  ad_ka_drag=0;
  ad_tau_wind=zeros(nx,2);
  ad_detady=zeros(nx,1);
  ad_d50=zeros(nx,1);
  ad_d90=zeros(nx,1);
  ad_pp=zeros(ns,1);
  ad_pp(i)=1;  % delta function
   

  % apply the full posterior covariance matrix to the adjoint vector, then
  % re-package the outputs for input to TL model
  phi=C2*[ad_h; ad_pp; ad_dgamma; ad_H0; ad_theta0; ad_ka_drag];
  ad_h=phi(1:nx);
  ad_params=struct;
  for j=1:length(fld)
    ad_params=setfield(ad_params,fld{j},phi(nx+j));
  end
  ad_dgamma = phi((nx+ns+1):(2*nx+ns));    % Change 20
  ad_H0=phi(nx+ns+nx+1);
  ad_theta0=phi(nx+ns+nx+2);
  ad_ka_drag=phi(nx+ns+nx+3);

  % tangent linear run, I*C*I'
  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp] = ...
      tl_hydroSedModel(ad_h,...
                       ad_H0,...
                       ad_theta0,...
                       0*ad_omega,...
                       ad_ka_drag,...
                       ad_dgamma,...
                       0*ad_tau_wind,...
                       0*ad_detady,...
                       0*ad_d50,...
                       0*ad_d90,...
                       ad_params,...
                       fcst);

  % extract the forecast bathy covariance for this gridpoint
  Chsp(:,nx+i)=[tl_hp; phi(nx+[1:ns])];

end  % loop over gridpoints
posterior.Chsp=Chsp;
