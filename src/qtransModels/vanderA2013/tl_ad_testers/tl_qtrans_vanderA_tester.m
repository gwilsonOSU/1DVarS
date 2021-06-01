%
% perturbation test of tl_qtrans_vanderA.m
%
clear

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
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
k=kabs*[cos(waves.theta(i)) sin(waves.theta(i))];
udelta=udelta_reniers2004(ubar,k,omega,h,...
                               Hrms,detady,windW,...
                               Dr,fv,d50);
ws=ws_brownLawler(d50);
param.n=1.2;
param.m=11;
param.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
param.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
param.fv=0.1;  % breaking-induced eddy viscosity calibration parameter, see
                % Reniers et al. (2004) Table 4.  Scalar, of order 0.1 (default)
[q,bkgd]=qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta,ws,param);

for niter=1:100

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_h          =bkgd.h          *frac_tl*myrand();
tl_d50        =bkgd.d50        *frac_tl*myrand();
tl_d90        =bkgd.d90        *frac_tl*myrand();
tl_h          =bkgd.h          *frac_tl*myrand();
tl_Hrms       =bkgd.Hrms       *frac_tl*myrand();
tl_kabs       =bkgd.kabs       *frac_tl*myrand();
tl_omega      =bkgd.omega      *frac_tl*myrand();
tl_udelta     =bkgd.udelta     *frac_tl*myrand();
tl_ws         =bkgd.ws         *frac_tl*myrand();
tl_param.n    =bkgd.param.n    *frac_tl*myrand();
tl_param.m    =bkgd.param.m    *frac_tl*myrand();
tl_param.xi   =bkgd.param.xi   *frac_tl*myrand();
tl_param.alpha=bkgd.param.alpha*frac_tl*myrand();
tl_param.fv   =bkgd.param.fv   *frac_tl*myrand();

% NL model with perturbation
param_prime.n    =param.n    +tl_param.n    ;
param_prime.m    =param.m    +tl_param.m    ;
param_prime.xi   =param.xi   +tl_param.xi   ;
param_prime.alpha=param.alpha+tl_param.alpha;
param_prime.fv   =param.fv   +tl_param.fv   ;
[q_prime,bkgd_prime] = ...
    qtrans_vanderA(d50+tl_d50,d90+tl_d90,h+tl_h,Hrms+tl_Hrms,...
                   kabs+tl_kabs,omega+tl_omega,udelta+tl_udelta,...
                   ws+tl_ws,param_prime);
tl_q_true(niter) = q_prime - q;

% TL model
tl_q(niter) = ...
    tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_Hrms,tl_kabs,...
                      tl_omega,tl_udelta,tl_ws,tl_param,bkgd);

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
