clear
close all

b0926 = matfile('bkgd0181_0926_1600EST.mat');

x_0926 = b0926.x;
d50_0926 = b0926.d50;
d90_0926 = b0926.d90;
h_0926 = b0926.h;
tanbeta_0926 = b0926.tanbeta;
Hrms_0926 = b0926.Hrms;
kabs_0926 = b0926.kabs;
omega_0926 = b0926.omega;
udelta_0926 = b0926.udelta_w;
delta_0926 = b0926.delta_bl;
ws_0926 = b0926.ws;
Aw_0926 = b0926.Aw;
Sw_0926 = b0926.Sw;
Uw_0926 = b0926.Uw;
Q_0926 = b0926.Q;
Q0_0926 = b0926.Q0;

%%
loc = 270;
pos_0926 = find(x_0926 == loc);
Q0 = Q0_0926(pos_0926);

param.alpha = 8.2; 
param.xi = 1.7;
param.m = 11;
param.n = 1.2;
param.nosusp = 1;

N = 1000;
n = 2500;

init = nan(N,n);
q_re = init;
q_mean = init;
q_bulk = [nan,nan];

t_est = N*n/2000*6/60;
prompt = ['Estimated runtime is ',num2str(t_est), ' minutes. Continue? [Y]'];
str = input(prompt,'s');

if str == 'Y'
    
    c = fix(clock);
    disp(['date: ',num2str(c(2)),'/',num2str(c(3))])
    disp(['time: ',num2str(c(4)),':',num2str(c(5))])
    tic
    
    d50 = d50_0926(pos_0926);
    d90 = d90_0926(pos_0926);
    h = h_0926(pos_0926) + b0926.tide;
    tanbeta = tanbeta_0926(pos_0926);
    kabs = kabs_0926(pos_0926);
    omega = omega_0926;
    udelta = udelta_0926(pos_0926,:);
    delta = delta_0926(pos_0926);
    ws = ws_0926(pos_0926);
    Hrms = Hrms_0926(pos_0926);

    k = kabs;
    B = Hrms/(sqrt(2));
    Hmo = 1.4*Hrms;
    [Aw,Sw,Uw,~] = Uwave_ruessink2012_params(Hmo,k,omega,h);
    q_bulk = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);

    parfor i = 1:N

        Omegatc = 0;

        for j = 1:n

            H = raylrnd(B);

            if H/h >= 0.78

                q_re(i,j) = nan;

            else

            Hrms = H;
            Hmo = H;

            [Aw,Sw,Uw,~] = Uwave_ruessink2012_params(Hmo,k,omega,h);
            [q_re(i,j),workspace_random,Omegatc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);

            end
        end

    end

    parfor i = 1:N
        for j = 1:n
            q_mean(i,j) = nanmean(q_re(i,1:j));
        end
    end

toc
    
else
    
end

s = std(q_mean);

save(['duck94_onshore_',num2str(loc),'.mat'],'q_bulk','q_mean','q_re','s')

figure()
histogram(abs(q_re))
set(gca,'xscale','log')
set(gca,'yscale','log')

figure()
txt = 'onshore';

for i = 1:N

    plot(q_mean(i,:),'LineWidth',2,'color',[0.8 0.8 0.8])
    hold on

end
plot(s,'b','LineWidth',4)

hold off

yline(Q0,'--k','LineWidth',2)
yline(q_bulk,'--b','LineWidth',2)
yline(0,'--k','LineWidth',2)
title('cumulative average transport')
subtitle(txt)
ylabel('q mean(m^2/s)')
xlabel('wave')
% set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([0 2500])
ylim([1E-7 1E-4])