clear 
close all

os_235 = matfile('duck94_onshore_235.mat');
b0926 = matfile('bkgd0181_0926_1600EST.mat');

q_bulk = os_235.q_bulk;
q_mean = os_235.q_mean;
q_re = os_235.q_re;
s = os_235.s;
%%
sz = size(q_mean);

omega = b0926.omega;
T = 2*pi/omega;
t = 0:T:T*(sz(2)-1);
t = t/3600;
t0 = 1.625E-6;
walk = t0*t.^-0.5;

figure('PaperPosition',[0 0 5 5])
hold on
plot(t,s,'-k','LineWidth',2)
plot(t, walk,'--k','LineWidth',2)
plot(t,mean(q_mean),'-r','LineWidth',2)
yline(q_bulk,'--b','LineWidth',2)

for i = 1:sz(1)
    plot(t,q_mean(i,:),'LineWidth',2,'color',[0.8 0.8 0.8])
    hold on
end
plot(t,s,'-k','LineWidth',2)
plot(t, walk,'--k','LineWidth',2)
plot(t,mean(q_mean),'-r','LineWidth',2)

hold off
yline(0,'--k','LineWidth',1)
title('x = 235m')
ylabel('q mean(m^2/s)')
xlabel('time (hrs)')
set(gca,'xscale','log')
set(gca,'FontSize',18)
xlim([0 6])
ylim([-1E-5 5E-5])
set(gca,'xtick',[1 3 6])
set(gca,'xticklabel',[1 3 6])
legend('std dev','t^{-1/2}','mean','q bulk')
print('cumulativemean235','-dpng','-r600')

%%
k = mean(kurtosis(q_re));
dim = [0.65 0.5 0.3 0.3];
str = ['kurtosis ',num2str(k)];

edgesp = 10.^(-6:0.1:-3);
[Np,edgesp] = histcounts(q_re,edgesp,'Normalization','countdensity');
edgesn = edgesp;
[Nn,edgesn] = histcounts(-q_re(q_re<0),edgesn,'Normalization','countdensity');
figure('PaperPosition',[0 0 5 4])
histx = linspace(2.5E-5,4E-4);

y03 = 0.5E-4;
y3 = y03*histx.^-3;
hold on
plot(histx,y3,'--k','LineWidth',2)
hold on
histogram('BinEdges',edgesn,'BinCounts',Nn,'FaceColor','r')
histogram('BinEdges',edgesp,'BinCounts',Np,'FaceColor','b')

set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',18)
set(gca,'xtick',[1E-6 1E-5 1E-4 1E-3])
set(gca,'ytick',[1E6 1E8 1E10 1E12])
hold on


hold off
title('x = 235m')
xlabel('q (m^2/s)')
ylabel('frequency')
legend('t^{-3}','offshore','onshore')
print('histogram235','-dpng','-r600')
