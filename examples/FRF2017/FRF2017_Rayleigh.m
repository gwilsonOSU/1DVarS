clear
close all
% comment
addpath(genpath('/Users/jack/Desktop/FRF2017'))
addpath(genpath('/Users/jack/Desktop/src'))

% load data
uE = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','currentEast');
p = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','pressure');
t = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','time');

% adjust t and u
t = t-t(1);
u = uE;
u_mean = mean(u);

% define h before detrending p
h = mean(p);

% filter u and p
fs = 2;
fcA = 1/15;  % wave band lower bound in Hz
fcB = 0.7;   % wave band upper bound in Hz
[bA,aA] = butter(5,fcA/(fs/2),'high');
[bB,aB] = butter(5,fcB/(fs/2),'high');
x = detrend(p);
xA = filtfilt(bA,aA,x);
pWB = xA - filtfilt(bB,aB,xA);
x = detrend(u);
xA = filtfilt(bA,aA,x);
uWB = xA - filtfilt(bB,aB,xA);

% find zero crossings (crude)
upzc = find(diff(sign(pWB))<0);

% calculate H
parfor i = 1:(length(upzc)-1)
    wavetime = t(upzc(i):upzc(i+1));
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1)); 
    
    r = corrcoef(wave_u,wave_p);
    r = abs(r(2,1));
    if r >= 0.9
        H(1,i) = max(wave_p)-min(wave_p);
        T(1,i) = wavetime(end)-wavetime(1);
    else
        H(1,i) = NaN;
        T(1,i) = NaN;
    end
end

% create pdf
H(find(isnan(H)))=[];
ray = fitdist(H.','Rayleigh');

% inputs to calculate q
T(find(isnan(T)))=[];
T_av = mean(T);
omega = 2*pi/T_av;
% c = sqrt(10*h);
% kabs = omega/c;
g = 9.81;
kabs = fzero( @(kabs) omega^2 - g*kabs*tanh(kabs*h), 0.1);
d50 = 0.00018;
d90 = 0.00024;
ws = ws_brownLawler(0.8*d50);
param.alpha = 8.2;
param.xi = 1.7;
param.m = 11;
param.n = 1.2;
udelta = [u_mean,0];

for i=1:2
    Hrms = random(ray);
    [q(1,i),k] = qtrans_vanderA_bulk(d50,d90,h,Hrms,kabs,omega,udelta,ws,param);
    if Hrms/h > 0.78
            breaking(1,i) = 1;
        else
            breaking(1,i) = 0;
        end
end

figure()
plot(q,'r')
title('net transport per wave for individual waves')
ylabel('q (m^2/s)')
xlabel('wave')
yline(0)

figure()
tiledlayout(3,1)
nexttile
histogram(abs(q),'FaceColor','r')
title('q, linear scale')
nexttile
histogram(abs(q),'FaceColor','r')
set(gca,'yscale','log')
title('q, log scale y axis')
nexttile
histogram(abs(q),'FaceColor','r')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('q, log log scale')

% HT = [H.' T.'];
% figure()
% hist3(HT,'FaceColor','r')
% xlabel('Wave height (m)')
% ylabel('Wave period (s)')
