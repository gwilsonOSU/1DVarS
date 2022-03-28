clear
close all

b0926 = matfile('bkgd0181_0926_1600EST.mat');

x_0926 = b0926.x;
h_0926 = b0926.h;
Hrms_0926 = b0926.Hrms;
Q_0926 = b0926.Q;
Q0_0926 = b0926.Q0;

yvars = [-h_0926, Q_0926, Hrms_0926];
titles = {'Bathymetry','Sediment transport rate (bulk wave)','Wave height'};
ylabels = {'h (m)','q (m^2/s)','Hrms (m)'};

figure('PaperPosition',[0 0 15 9])
for i = 1:3
    subplot(3,1,i)
    plot(x_0926,yvars(:,i)','LineWidth',3)
    xline(195,'--k',{'a'},'LineWidth',2)
    xline(235,'--k',{'b'},'LineWidth',2)
    xline(270,'--k',{'c'},'LineWidth',2)
    if i == 1
        ylim([-4 0])
    else
    end
    xlim([125 305])
    title(titles(i))
    ylabel(ylabels(i))
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',24)
    xlabel('Cross shore distance (m)')
end
set(gca,'xtick',[125,160,195,235,270,305])
set(gca,'xticklabel',[125,160,195,235,270,305])
print('hhq','-dpng','-r600')