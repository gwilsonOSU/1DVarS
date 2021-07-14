clear

dir = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/c/Trial2/';


load(strcat(dir,'Time.mat'));
TimeTrial = Time;


load(strcat(dir,'RMSE.mat'));
RMSETrial = RMSE;

load(strcat(dir,'RMSEp.mat'));
RMSEpTrial = RMSEp;

%Plotting RMSE first:
figure('Units', 'inches','Position', [0,0,15,7])
plot(TimeTrial, RMSETrial, '-og')
% legend ('Test a','Test b','Test c','FontSize',14,'Interpreter','Latex','Location','northwest')
xlabel ('Time (hrs) from 03-Oct-1994 07:00:16 EST','FontSize',14,'Interpreter','Latex')
ylabel ('RMSE (${m^2}$)','FontSize',14,'Interpreter','Latex')
title('RMS Error in assimilated bathymetry - Test C','FontSize',14,'Interpreter','Latex')
grid on
print ('-dpng', '-r300', strcat(dir, 'RMSE','.png'))
close
 
%Time 1: 01-Sep-1994 16:00:16  datenum: 728538.666851852
% Time 3: Time + 03-Oct-1994 07:00:16  (datenum: 728570.291851852)
%Time 2:  24-Sep-1994 16:00:16 datenum: 728561.666851852



%Plotting RMSEp
figure('Units', 'inches','Position', [0,0,15,7])
plot(TimeTrial, RMSEpTrial, '-og')
% legend ('Test a','Test b','Test c','FontSize',14,'Interpreter','Latex','Location','northwest')
xlabel ('Time (hrs) from 03-Oct-1994 07:00:16 EST','FontSize',14,'Interpreter','Latex')
ylabel ('RMSE (${m^2}$)','FontSize',14,'Interpreter','Latex')
title('RMS Error in forecast bathymetry - Test C','FontSize',14,'Interpreter','Latex')
grid on
print ('-dpng', '-r300', strcat(dir, 'RMSEp','.png'))
close









