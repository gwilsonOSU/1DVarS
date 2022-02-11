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

figure()
for i = 1:3
    subplot(3,1,i)
    plot(x_0926,yvars(:,i)','LineWidth',2)
    xline(195,'--k',{'a'})
    xline(235,'--k',{'b'})
    xline(270,'--k',{'c'})
    xlim([100 400])
    title(titles(i))
    ylabel(ylabels(i))
    set(gca,'xtick',[90,125,160,195,235,270,305,340,375])
    set(gca,'xticklabel',[90,125,160,195,235,270,305,340,375])
    xlabel('Cross shore distance (m)')
end
