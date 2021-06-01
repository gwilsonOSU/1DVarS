%
% perturbation test of tl_Uwave_ruessink2012.m
%
clear

load ../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=50;
bkgd.Hmo=waves.H(i)*1.4;
bkgd.k=waves.k(i);
bkgd.omega=waves.sigma;
bkgd.h=waves.h(i);
phs=linspace(0,2*pi,1000);

[u,bkgd_uwave] = ...
    Uwave_ruessink2012(phs,...
                       bkgd.Hmo,...
                       bkgd.k,...
                       bkgd.omega,...
                       bkgd.h);

for niter=1:100

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_Hmo  =bkgd.Hmo  -bkgd.Hmo  *(frac_tl*myrand()+1);
tl_k    =bkgd.k    -bkgd.k    *(frac_tl*myrand()+1);
tl_omega=bkgd.omega-bkgd.omega*(frac_tl*myrand()+1);
tl_h    =bkgd.h    -bkgd.h    *(frac_tl*myrand()+1);

% compare TL perturbation to actual
u_prime = ...
    Uwave_ruessink2012(phs,...
                       bkgd.Hmo+tl_Hmo,...
                       bkgd.k+tl_k,...
                       bkgd.omega+tl_omega,...
                       bkgd.h+tl_h);
tl_u(:,niter) = tl_Uwave_ruessink2012(tl_Hmo,tl_k,tl_omega,tl_h,phs,bkgd_uwave);

tl_u_true(:,niter) = u_prime - u;

end

% Compare.  At first I was plotting percentage errors, but this
% turned out to be flawed as I was seeing a lot of division by small
% numbers.  Since then I changed to a scatterplot comparison, which is a
% much better indicator of validity
clf
itest=500;
plot(tl_u_true(itest,:),tl_u(itest,:),'.')
hold on
axis equal tight
ax=axis;
plot(ax([1 2]),ax([1 2]),'k--')
hold off
title(['u comparison (should be 1-1)'])
xlabel('true [m/s]')
ylabel('predicted [m/s]')

% clf
% plot(tl_Uw_true,tl_Uw,'.')
% hold on
% axis equal tight
% ax=axis;
% plot(ax([1 2]),ax([1 2]),'k--')
% hold off
