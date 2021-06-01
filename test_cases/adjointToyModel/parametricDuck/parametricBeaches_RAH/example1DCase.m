% example1DCase, compute parametric bar predictions for a test case at Duck.  
% Compare measured an example CRAB survey survey.  This is for 1D transects
% only.
% Holman, May, 2016

clear
load('testCase1D')

% pick an example profile from the CRAB
i = 81;                     % arbitrary example, chose any case you want.
figure(1); clf
plot(survey.x,survey.Z(i,:))
grid on; xlabel('x (m)'); ylabel('z (m)')

% manually chose a shoreline, then the sand bar locations
disp('Digitize the shoreline location (h = 0)')
[xs,~] = ginput(1);
betaShore = 0.1;            % average shoreline beach slope (from literature)
disp('Now digitize the sand bar crest location');
[xb,~] = ginput(1);

xOff = 700; hOff = 7.5;     % chose a deep point from other information
betaOff = 0.0095;           % bathy slope at deep point
hSea = 4.5;                 % from Ruessink paper for Duck
h = make1DBeachEngine(survey.x,xs,betaShore,xb,xOff,hOff,betaOff,hSea);

figure(1); hold on
plot(survey.x,-h,'r')
legend('survey', 'parametric')
