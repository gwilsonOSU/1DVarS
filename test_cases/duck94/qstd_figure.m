close all
clear

qstd = matfile('Qstd0181_0926_1600EST.mat');

Qstd = qstd.Qstd;
x = qstd.x;

figure('PaperPosition',[0 0 10 4])
plot(x,Qstd,'LineWidth',3)
xlim([125 305])
set(gca,'xtick',[125,160,195,235,270,305])
set(gca,'xticklabel',[125,160,195,235,270,305])
set(gca,'FontSize',18)
xline(195,'--k','LineWidth',2)
xline(235,'--k','LineWidth',2)
xline(270,'--k','LineWidth',2)
xlabel('x (m)')
ylabel('std dev of Q (m^2/s)')

print('qstdx','-dpng','-r600')