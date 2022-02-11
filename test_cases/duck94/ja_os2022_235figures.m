clear 
close all

os_235 = matfile('duck94_onshore_235.mat');
b0926 = matfile('bkgd0181_0926_1600EST.mat');

q_bulk = os_235.q_bulk;
q_mean = os_235.q_mean;
q_re = os_235.q_re;
s = os_235.s;

sz = size(q_mean);

omega = b0926.omega;
T = 2*pi/omega;
t = 0:T:T*(sz(2)-1);
t = t/3600;
t0 = 1.625E-6;
walk = t0*t.^-0.5;

figure()
hold on
plot(t,s,'-k','LineWidth',2)
plot(t, walk,'--k','LineWidth',2)
plot(t,mean(q_mean),'-r','LineWidth',2)

for i = 1:sz(1)
    plot(t,q_mean(i,:),'LineWidth',2,'color',[0.8 0.8 0.8])
    hold on
end
plot(t,s,'-k','LineWidth',2)
plot(t, walk,'--k','LineWidth',2)
plot(t,mean(q_mean),'-r','LineWidth',2)


hold off
yline(q_bulk,'--b',{'q bulk'},'LineWidth',2)
yline(0,'--k','LineWidth',1)
title('Cumulative average transport - bar crest (x = 235m)')
ylabel('q mean(m^2/s)')
xlabel('time (hrs)')
set(gca,'xscale','log')
% set(gca,'yscale','log')
xlim([0 6])
ylim([-1E-5 5E-5])
set(gca,'xtick',[0.5 1 3 6])
set(gca,'xticklabel',[0.5 1 3 6])
legend('standard deviation','t^{-1/2}','mean')


%%
k = mean(kurtosis(q_re));
dim = [0.65 0.5 0.3 0.3];
str = ['kurtosis ',num2str(k)];

edges = 10.^(-6:0.1:-3);
[N,edges] = histcounts(abs(q_re),edges,'Normalization','countdensity');

figure()
histx = linspace(2.5E-5,4E-4);
y02 = 0.2E1;
y2 = y02*histx.^-2;

plot(histx,y2,'--k','LineWidth',2)
hold on
histogram('BinEdges',edges,'BinCounts',N,'FaceColor','b')
set(gca,'yscale','log')
set(gca,'xscale','log')
hold on
annotation('textbox',dim,'String',str,'FitBoxToText','on')
hold off
title('Distribution of magnitude of transport - bar crest (x = 235m)')
xlabel('q (m^2/s)')
ylabel('frequency')
legend('slope -2')
% figure()
% histogram(abs(q_re))
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% hold on
% 
% histx = linspace(4E-5,1E-3);
% y02 = 4E-7;
% y2 = y02*histx.^-2;
% y03 = 1E-11;
% y3 = y03*histx.^-3;
% 
% plot(histx,y2,'--b','LineWidth',2)
% plot(histx,y3,'--b','LineWidth',2)
