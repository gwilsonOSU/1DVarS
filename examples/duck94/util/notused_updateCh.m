
% OLD: this code used to be in assim_1dh.m, but became mostly obsolete.
% Moved to here to clean things up.

% calculate covariance of the forecast bathymetry, and its covariance with
% sediment transport parameters.  This will facilitate parameter corrections
% in future analysis cycles
if(~doCovUpdate)  % skip this step
  disp('not calculating forecast covariance Chp')
  posterior.Chp = prior.Ch;
  return;
end
disp('propagating forecast covariance')
Chp=zeros(nx);
parfor i=1:nx
  % disp(['  gridpoint ' num2str(i) ' of ' num2str(nx)])
  if(floor(i/nx*10)>floor((i-1)/nx*10))
    disp([num2str(floor(i/nx*10)*10) '%'])
  end

  % adjoint model acting on identity matrix
  zz=zeros(nx,nsubsteps);
  comb=zz;
  comb(i,end)=1;  % data functional (delta-fn, aka identity matrix)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,ad_dgamma,...
   ad_tau_wind,ad_detady,ad_d50,ad_d90,ad_params] = ...
      ad_hydroSedModel(zz,zz,zz,zz,zz,comb,bkgd);

  % apply covariance to model inputs
  phi=C2*[ad_h; ad_dgamma; ad_H0; ad_theta0; ad_ka_drag];
  ad_h=phi(1:nx);
  ad_dgamma =phi(nx+[1:nx]);
  ad_H0     =phi(2*nx+1);
  ad_theta0 =phi(2*nx+2);
  ad_ka_drag=phi(2*nx+3);

  % apply covariance to sed params
  ad_pp=zeros(ns,1);  % init
  for j=1:ns
    ad_pp(j)=getfield(ad_params,sedParamList{j});  
  end
  ad_pp = Cs*ad_pp;
  for j=1:ns
    ad_params=setfield(ad_params,sedParamList{j},ad_pp(j));
  end

  % tangent linear run, I*C*I'
  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp] = ...
    tl_hydroSedModel(ad_h,...
                     ad_H0,...
                     ad_theta0,...
                     0*ad_omega,...
                     ad_ka_drag,...
                     0*ad_tau_wind,...
                     0*ad_detady,...
                     ad_dgamma,...
                     0*ad_d50,...
                     0*ad_d90,...
                     ad_params,...
                     bkgd);

  % extract the forecast bathy covariance for this gridpoint
  Chp(:,i)=tl_hp(:,end);

end  % loop over gridpoints
posterior.Chp=Chp;
