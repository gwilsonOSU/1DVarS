%Script to plot errors in assimilation and forecast:
clear

dira = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/a/Trial2/';
dirb = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/b/Trial2/';
dirc = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/c/Trial2/';

% dir3 = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial3/';
% dir7 = '/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Trials/Trial7/';


load(strcat(dira,'Time.mat'));
TimeTriala = Time;
load(strcat(dirb,'Time.mat'));
TimeTrialb = Time;
load(strcat(dirc,'Time.mat'));
TimeTrialc = Time;







load(strcat(dira,'RMSE.mat'));
RMSETriala = RMSE;
load(strcat(dirb,'RMSE.mat'));
RMSETrialb = RMSE;
load(strcat(dirc,'RMSE.mat'));
RMSETrialc = RMSE;







load(strcat(dira,'RMSEp.mat'));
RMSEpTriala = RMSEp;
load(strcat(dirb,'RMSEp.mat'));
RMSEpTrialb = RMSEp;
load(strcat(dirc,'RMSEp.mat'));
RMSEpTrialc = RMSEp;



%Plotting RMSE first:
figure('Units', 'inches','Position', [0,0,15,7])
plot(TimeTriala, RMSETriala, '-og')
hold on
plot(TimeTrialb, RMSETrialb, '-ob')
hold on
plot(TimeTrialc, RMSETrialc, '-or')
hold on
legend ('Test a','Test b','Test c','FontSize',14,'Interpreter','Latex','Location','northwest')
xlabel ('Time (hrs)','FontSize',14,'Interpreter','Latex')
ylabel ('RMSE (${m^2}$)','FontSize',14,'Interpreter','Latex')
title('RMS Error in assimilated bathymetry','FontSize',14,'Interpreter','Latex')

% print ('-dpng', '-r300', strcat('/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Plots/RMSEforTrials','.png'))

% Time 3: Time + 




%Plotting RMSE forecast:
figure('Units', 'inches','Position', [0,0,15,7])
plot(TimeTriala, RMSEpTriala, '-og')
hold on
plot(TimeTrialb, RMSEpTrialb, '-ob')
hold on
plot(TimeTrialc, RMSEpTrialc, '-or')
legend ('Test a','Test b','Test c','FontSize',14,'Interpreter','Latex','Location','northwest')
xlabel ('Time (hrs)','FontSize',14,'Interpreter','Latex')
ylabel ('RMSE Forecast (${m^2}$)','FontSize',14,'Interpreter','Latex')
title('RMS Error in Forecasts w/ Sediment Transportation','FontSize',14,'Interpreter','Latex')
% print ('-dpng', '-r300', strcat('/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Plots/RMSE_p_forTrials','.png'))




































% Plotting RMSE for different trials for period b:
% 
% load(strcat(dir1,'Time.mat'));
% TimeTrial1 = Time;
% load(strcat(dir2,'Time.mat'));
% TimeTrial2 = Time;
% load(strcat(dir4,'Time.mat'));
% TimeTrial4 = Time;
% load(strcat(dir3,'Time.mat'));
% TimeTrial3 = Time;
% load(strcat(dir7,'Time.mat'));
% TimeTrial7 = Time;
% 
% 
% 
% 
% 
% 
% load(strcat(dira,'RMSE.mat'));
% RMSETrial1 = RMSE;
% load(strcat(dirb,'RMSE.mat'));
% RMSETrial2 = RMSE;
% load(strcat(dirc,'RMSE.mat'));
% RMSETrial4 = RMSE;
% load(strcat(dir3,'RMSE.mat'));
% RMSETrial3 = RMSE;
% load(strcat(dir7,'RMSE.mat'));
% RMSETrial7 = RMSE;
% 
% 
% 
% 
% 
% 
% load(strcat(dira,'RMSEp.mat'));
% RMSEpTrial1 = RMSEp;
% load(strcat(dirb,'RMSEp.mat'));
% RMSEpTrial2 = RMSEp;
% load(strcat(dirc,'RMSEp.mat'));
% RMSEpTrial4 = RMSEp;
% load(strcat(dir3,'RMSEp.mat'));
% RMSEpTrial3 = RMSEp;
% load(strcat(dir7,'RMSEp.mat'));
% RMSEpTrial7 = RMSEp;
% 
% %Plotting RMSE first:
% figure('Units', 'inches','Position', [0,0,15,7])
% plot(TimeTrial1, RMSETrial1, '-or')
% hold on
% plot(TimeTrial2, RMSETrial2, '-ob')
% hold on
% plot(TimeTrial3, RMSETrial3, '-og')
% hold on
% plot(TimeTrial4, RMSETrial4, '-oc')
% hold on
% plot(TimeTrial7, RMSETrial7, '-ok')
% legend ('Trial1 : No CoV propagation ','Trial2: Forecast CoV propagation','Trial3 cgamma = 0','Trial4 params.std = 0','Trial7 max Ch = 0.5','FontSize',14,'Interpreter','Latex','Location','northwest')
% xlabel ('Time (hrs)','FontSize',14,'Interpreter','Latex')
% ylabel ('RMSE (${m^2}$)','FontSize',14,'Interpreter','Latex')
% title('RMS Error in assimilated bathymetry','FontSize',14,'Interpreter','Latex')
% 
% print ('-dpng', '-r300', strcat('/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Plots/RMSEforTrials','.png'))
% 
% 
% 
% 
% 
% 
% %Plotting RMSE forecast:
% figure('Units', 'inches','Position', [0,0,15,7])
% plot(TimeTrial1, RMSEpTrial1, '-or')
% hold on
% plot(TimeTrial2, RMSEpTrial2, '-ob')
% hold on
% plot(TimeTrial3, RMSEpTrial3, '-og')
% hold on
% plot(TimeTrial4, RMSEpTrial4, '-oc')
% hold on
% plot(TimeTrial7, RMSEpTrial7, '-ok')
% legend ('Trial1 : No CoV propagation ','Trial2: Forecast CoV propagation','Trial3 cgamma = 0','Trial4 params.std = 0','Trial7 max Ch = 0.5','FontSize',14,'Interpreter','Latex','Location','northwest')
% xlabel ('Time (hrs)','FontSize',14,'Interpreter','Latex')
% ylabel ('RMSE Forecast (${m^2}$)','FontSize',14,'Interpreter','Latex')
% title('RMS Error in Forecasts w/ Sediment Transportation','FontSize',14,'Interpreter','Latex')
% print ('-dpng', '-r300', strcat('/home/server/student/homes/asalim/Desktop/DataAssimilation/1DVarS-main/test_cases/duck94/Plots/RMSE_p_forTrials','.png'))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Plotting Ch:
% % figure('Units', 'inches','Position', [0,0,15,7])
% % plot(TimeTrial1, ChTrial1, '-og')
% % hold on
% % plot(TimeTrial2, ChTrial2, '-ob')
% % hold on
% % plot(TimeTrial3, ChTrial3, '-or')
% % legend ('Trial1 : No CoV propagation ','Trial2: Forecast CoV propagation','Trial 3: Forecast CoV propagation w/ CGamma = 0','FontSize',14,'Interpreter','Latex','Location','northwest')
% % xlabel ('Time (hrs)','FontSize',14,'Interpreter','Latex')
% % ylabel ('Ch E(i,i)(${m^2}$)','FontSize',14,'Interpreter','Latex')
% % title('Ch variation ','FontSize',14,'Interpreter','Latex')
% 
% 
% 
% 
% 
% 
% 









