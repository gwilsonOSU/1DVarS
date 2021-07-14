%Plotting Ch Maxs:
clear


directory = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6/Ch/';
files = dir(strcat(directory,'*.mat'));

n = 1:1:length(files);
ChMax = zeros(length(n),1);

% XFRF = 100:5:900;
% XFRF = flip(XFRF);
% directoryPlot = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6/PcolorPlots/';

for i = 1:length(files)
    mat = load(strcat(directory,'Ch',num2str(i)));
    Ch = mat.Ch;
    ChMax(i,1) = max(max(Ch));
end

figure('Units', 'inches','Position', [0,0,8,4])
plot(n, ChMax, '-*b')
xlabel ('n','FontSize',14,'Interpreter','Latex')
ylabel ('Ch max (m^2)','FontSize',14,'Interpreter','Latex')
title('Ch max variation ','FontSize',14,'Interpreter','Latex')
grid on
print ('-dpng', '-r300', strcat('/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6/','ChMAx ','.png'))



%     s = pcolor(XFRF,XFRF,Ch);
%     title(strcat('Ch for n = ',num2str(i)),'FontSize',14,'Interpreter','Latex')
%     s.FaceColor = 'interp';
%     shading flat
%     c = colorbar;
%     % caxis([-e-3 e-3])
%     dc = (max(max(Ch)) - min(min(Ch)))/10;
%     c.Ticks = [min(min(Ch)): dc: max(max(Ch))];
%     colormap jet     
%     print ('-dpng', '-r300', strcat(directoryPlot, 'Ch- ','n= ', num2str(i),'.png'))
%     close





















