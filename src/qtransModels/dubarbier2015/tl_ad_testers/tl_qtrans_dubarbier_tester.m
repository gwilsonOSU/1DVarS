%
% perturbation test of tl_qtrans_dubarbier.m
%
clear

load ../../../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;

% choose a gridpoint and make a realistic background state
i=25;
bkgd.h=waves.h(i);
bkgd.Hrms=waves.H(i);
bkgd.ubar(1)=-(waves.E(i)+2*waves.Er(i))./(waves.c(i).*waves.h(i));
bkgd.ubar(2)=waves.v(i);
bkgd.k=waves.k(i);
bkgd.omega=waves.sigma;
bkgd.detady=.001;
bkgd.windW =[1 1]*.5;
bkgd.Dr=waves.eps_r(i)*1030;
bkgd.fv=1;
bkgd.d50=200e-6;
bkgd.k=bkgd.k*[cos(waves.theta(i)) sin(waves.theta(i))];
bkgd.udelta=udelta_reniers2004(bkgd.ubar,bkgd.k,bkgd.omega,bkgd.h,...
                               bkgd.Hrms,bkgd.detady,bkgd.windW,...
                               bkgd.Dr,bkgd.fv,bkgd.d50);
bkgd.ws=ws_brownLawler(bkgd.d50);
bkgd.Cw=0.00483;
bkgd.Cc=0.02002;
bkgd.Cf=0.01173;
bkgd.Ka=0.631e-4;
bkgd.tanbeta=0.01;
if(length(i)==1)
  bkgd.x=bkgd.tanbeta;  % wonky code hack
end

for niter=1:100

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_tanbeta=bkgd.tanbeta-bkgd.tanbeta*(frac_tl*myrand()+1);
tl_h      =bkgd.h      -bkgd.h      *(frac_tl*myrand()+1);
tl_Hrms   =bkgd.Hrms   -bkgd.Hrms   *(frac_tl*myrand()+1);
tl_k      =bkgd.k      -bkgd.k      *(frac_tl*myrand()+1);
tl_omega  =bkgd.omega  -bkgd.omega  *(frac_tl*myrand()+1);
tl_udelta =bkgd.udelta -bkgd.udelta *(frac_tl*myrand()+1);
tl_ws     =bkgd.ws     -bkgd.ws     *(frac_tl*myrand()+1);
tl_Cw     =bkgd.Cw     -bkgd.Cw     *(frac_tl*myrand()+1);
tl_Cc     =bkgd.Cc     -bkgd.Cc     *(frac_tl*myrand()+1);
tl_Cf     =bkgd.Cf     -bkgd.Cf     *(frac_tl*myrand()+1);
tl_Ka     =bkgd.Ka     -bkgd.Ka     *(frac_tl*myrand()+1);

% test code: cancel some of them
% tl_tanbeta=0;
% tl_h      =0*tl_h      ;
% tl_Hrms   =0*tl_Hrms   ;
% tl_k      =0*tl_k      ;
% tl_udelta =0*tl_udelta ;
% tl_ws     =0*tl_ws     ;
% tl_Cw     =0*tl_Cw     ;
% tl_Cc     =0*tl_Cc     ;
% tl_Cf     =0*tl_Cf     ;
% tl_Ka     =0*tl_Ka     ;

% compare TL perturbation to actual
[q,qb,qs,qa,bkgd_qt]=qtrans_dubarbier(bkgd.tanbeta,...
                                      bkgd.h      ,...
                                      bkgd.Hrms   ,...
                                      bkgd.k      ,...
                                      bkgd.omega  ,...
                                      bkgd.udelta ,...
                                      bkgd.ws     ,...
                                      bkgd.Cw     ,...
                                      bkgd.Cc     ,...
                                      bkgd.Cf     ,...
                                      bkgd.Ka     );
q_prime = ...
    qtrans_dubarbier(bkgd.tanbeta+tl_tanbeta,...
                     bkgd.h      +tl_h      ,...
                     bkgd.Hrms   +tl_Hrms   ,...
                     bkgd.k      +tl_k      ,...
                     bkgd.omega  +tl_omega  ,...
                     bkgd.udelta +tl_udelta ,...
                     bkgd.ws     +tl_ws     ,...
                     bkgd.Cw     +tl_Cw     ,...
                     bkgd.Cc     +tl_Cc     ,...
                     bkgd.Cf     +tl_Cf     ,...
                     bkgd.Ka     +tl_Ka     );
tl_q(niter) = ...
    tl_qtrans_dubarbier(tl_tanbeta,...
                        tl_h,...
                        tl_Hrms,...
                        tl_k,...
                        tl_omega,...
                        tl_udelta,...
                        tl_ws,...
                        tl_Cw,...
                        tl_Cc,...
                        tl_Cf,...
                        tl_Ka,...
                        bkgd_qt);
tl_q_true(niter) = q_prime - q;

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
