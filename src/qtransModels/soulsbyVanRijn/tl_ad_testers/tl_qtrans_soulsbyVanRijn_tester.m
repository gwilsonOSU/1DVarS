%
% perturbation test of tl_qtrans_vanderA.m
%
clear

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
x=waves.x(i);
h=waves.h(i);
Hrms=waves.H(i);
ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
ubar(2)=waves.v(i);
kabs=waves.k(i);
omega=waves.sigma;
detady=.001;
windW =[1 1]*.5;
Dr=waves.eps_r(i)*1030;
fv=1;
d50=180e-6;
d90=400e-6;
ws=ws_brownLawler(d50);
param.alphab=1.6;
param.facua =1.0;
tanbeta=.01;  % scalar
theta=waves.theta(i);

[q,bkgd] = qtrans_soulsbyVanRijn(x,d50,d90,h,tanbeta,Hrms,kabs,...
                                 omega,theta,ubar,Dr,param);

for niter=1:100

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_d50    =bkgd.d50    *frac_tl*myrand();
tl_d90    =bkgd.d90    *frac_tl*myrand();
tl_h      =bkgd.h      *frac_tl*myrand();
tl_tanbeta=bkgd.tanbeta*frac_tl*myrand();
tl_Hrms   =bkgd.Hrms   *frac_tl*myrand();
tl_kabs   =bkgd.kabs   *frac_tl*myrand();
tl_omega  =bkgd.omega  *frac_tl*myrand();
tl_theta  =bkgd.theta  *frac_tl*myrand();
tl_ubar   =bkgd.ubar   *frac_tl*myrand();
tl_Dr     =bkgd.Dr     *frac_tl*myrand();
tl_param.alphab=param.alphab*frac_tl*myrand();
tl_param.facua =param.facua *frac_tl*myrand();

% NL model with perturbation
param_prime.alphab=param.alphab+tl_param.alphab;
param_prime.facua =param.facua +tl_param.facua ;
q_prime = qtrans_soulsbyVanRijn(x      ,...
                                d50    +tl_d50    ,...
                                d90    +tl_d90    ,...
                                h      +tl_h      ,...
                                tanbeta+tl_tanbeta,...
                                Hrms   +tl_Hrms   ,...
                                kabs   +tl_kabs   ,...
                                omega  +tl_omega  ,...
                                theta  +tl_theta  ,...
                                ubar   +tl_ubar   ,...
                                Dr     +tl_Dr     ,...
                                param_prime);
tl_q_true(niter) = q_prime - q;

% TL model
tl_q(niter) = tl_qtrans_soulsbyVanRijn(tl_d50,tl_d90,tl_h,tl_tanbeta,...
                                       tl_Hrms,tl_kabs,tl_omega,tl_theta,tl_ubar,...
                                       tl_Dr,tl_param,bkgd);

end

% show all stats.  At first I was plotting percentage errors, but this
% turned out to be flawed as I was seeing a lot of division by small
% numbers.  Since then I changed to a scatterplot comparison, which is a
% much better indicator of validity
clf
plot(tl_q_true,tl_q,'.')
hold on
axis equal tight
ax=axis;
plot(ax([1 2]),ax([1 2]),'k--')
hold off
xlabel('true')
ylabel('predicted')
