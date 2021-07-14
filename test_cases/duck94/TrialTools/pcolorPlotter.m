clear
% cd /home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6

directory = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6/Ch/';
files = dir(strcat(directory,'*.mat'));
XFRF = 100:5:900;
XFRF = flip(XFRF);
directoryPlot = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial6/PcolorPlots/';


parfor i = 1:length(files)
    mat = load(strcat(directory,files(i).name));
    Ch = mat.Ch;
    s = pcolor(XFRF,XFRF,Ch);
    title(strcat('Ch for n = ',num2str(i)),'FontSize',14,'Interpreter','Latex')
    s.FaceColor = 'interp';
    shading flat
    c = colorbar;
    % caxis([-e-3 e-3])
    dc = (max(max(Ch)) - min(min(Ch)))/10;
    c.Ticks = [min(min(Ch)): dc: max(max(Ch))];
    colormap jet     
    print ('-dpng', '-r300', strcat(directoryPlot, 'Ch- ','n= ', num2str(i),'.png'))
    close
end



























