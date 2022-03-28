%% file loading stuff
clear 
close all

addpath(genpath('/Users/jack/Documents/Github/1DVarS'))
% addpath(genpath('/home/shusin8/users/wilsongr/Duck94_SteveH/'))
% addpath(genpath('/home/server/student/homes/aldricja/1DVarS'))

m = matfile('0926.mat');
dat = m.dat;
pos = m.pos;

%% define variables
fs = 2;
time = 6; % time in 8hr blocks 
loc = 10;

u_pos = pos(time,2,loc,1);
p_pos = pos(time,1,loc,1);
u = dat(time,2,loc,:);
u = u(:);
p = dat(time,1,loc,:);
p = p(:);
t = 0:1/fs:(length(u)-1)/2;
t = t';

%% check plot before filtering
close all

p_plot = detrend(p);

figure()
plot(t,u,t,p_plot)
legend('velocity','pressure')
xlim([0 100])
title('sanity check')
%% filter u and p
% define moving mean for u and p
filter_time = 512;
u_mean = movmean(u,filter_time);
p_mean = movmean(p, filter_time);

% enter wave band filter parameters
fcA = 1/20;  % wave band lower bound in Hz
fcB = 0.5;   % wave band upper bound in Hz
[bA,aA] = butter(5,fcA/(fs/2),'high');
[bB,aB] = butter(5,fcB/(fs/2),'high');
[bC,aC] = butter(5,fcB/(fs/2),'low');
% detrend and filter u
% wave band filter u 
x = detrend(u);
xA = filtfilt(bA,aA,x);
uWB = xA - filtfilt(bB,aB,xA);

% low pass filter u
x = detrend(u);
xA = filtfilt(bC,aC,x);
uLP = xA;

% detrend and filter p
x = detrend(p);
xA = filtfilt(bA,aA,x);
pWB = xA - filtfilt(bB,aB,xA);

%% check correlation before calculating q
% Chop into waves roughly based on p zero crossings
upzc = find(diff(sign(pWB))>0);
init = NaN(1,length(upzc)-1);
correlated = init;
% loop to check correlation
for i = 1:(length(upzc)-1)
    % grabbing rough t, u, and p for each wave
    wave_t = t(upzc(i):upzc(i+1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% possible bug!! should be t(upzc(i)+1) ???
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1)); 
    
    % checking correlation between p and u and calculating other inputs
    r = corrcoef(wave_u,wave_p);
    r = abs(r(2,1));
    if r >= 0.9
        correlated(1,i) = 1;
    else
        correlated(1,i) = 0;
    end
end

%% calculate q using pressure whenever possible
% Chop into waves roughly based on p zero crossings
upzc = find(diff(sign(pWB))>0);

% input parameters relevant to all waves
d50 = 0.00018;
d90 = 0.00024;
wsc = ws_brownLawler(0.8*d50);
wst = wsc;
ws = wsc;

param.alpha = 8.2;
param.xi = 1.7;
param.m = 11;
param.n = 1.2;
g = 9.81;

% initializing vectors for other inputs
q_vec = init;
H_vec = init;
T_vec = init;
breaking = init;
q_mean = init;
weird_wave = init;
q_param_vec = init;
q_mean_param = init;

cell_init = cell(1,length(upzc)-1);
u_param = cell_init;
t_param = cell_init;
udelta_cell = cell_init;

Omegatc = 0;
% loop to calculate other inputs and q
for i = 1:(length(upzc)-1)
    % grabbing rough t, u, and p for each wave
    wave_t = t(upzc(i):upzc(i+1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% possible bug!! should be t(upzc(i)+1) ???
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1)); 
    
    % checking correlation between p and u and calculating other inputs
    if correlated(1,i) == 1
        % udelta as mean of u at rough start of wave 
        udelta = [mean(wave_u),0];
        wave_u = wave_u - udelta(1);
        
        % uhat, uhatc, and uhatt based on definitions
        uhat = sqrt(2*mean(wave_u.^2));        
        uhatc = max(wave_u);        
        uhatt = abs(min(wave_u));
        
        % interpolated t for beginning of wave
        x10 = t(upzc(i));
        x20 = t(upzc(i)+1);
        y10 = pWB(upzc(i));
        y20 = pWB(upzc(i)+1);
        m0 = (y20-y10)/(x20-x10);
        t00 = x10 - y10/m0; 
        
        % interpolating t for end of wave
        x12 = t(upzc(i+1));
        x22 = t(upzc(i+1)+1);
        y12 = pWB(upzc(i+1));
        y22 = pWB(upzc(i+1)+1);
        m2 = (y22-y12)/(x22-x12);
        t02 = x12 - y12/m2;
        
        % calculating T
        T = t02 - t00;
        T_vec(1,i) = T;
        
        % interpolating t for mid point of wave
        midzc = find(diff(sign(wave_p))<0);
        x11 = wave_t(midzc);
        x21 = wave_t(midzc+1);
        y11 = wave_p(midzc);
        y21 = wave_p(midzc+1);
        m1 = (y21-y11)/(x21-x11);
        t01 = x11 - y11/m1;
        
        % calculating Tc and Tt
        Tc = t01 - t00;
        Tt = t02 - t01;
        
        % calculating Ttu and Tcu, with precautions to keep Ttu positive
        zcross = floor(t01);
        trough = find(wave_t>=zcross);
        [maxcrest,maxcrestidx] = max(wave_p);
        [mintrough,mintroughidx] = min(wave_p(trough));
        
        Ttu = wave_t(mintroughidx+length(wave_t)-length(trough)) - t01;        
        Tcu = wave_t(maxcrestidx) - t00;
        
        % dispersion solver for c
        h = p_mean(upzc(i));        
        omega = 2*pi/T;
        k = fzero( @(k) omega^2 - g*k*tanh(k*h), 0.1);
        kabs = k;
        c = omega/k;
        
        % calculate q
        [q_vec(1,i),~,Omegatc] = qtrans_vanderA_onewave_ja(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param,Omegatc);
        q_mean(1,i) = nanmean(q_vec);
        % calculate wave height
        H = max(wave_p) - min(wave_p);
        H_vec(1,i) = H;
        Hrms = H;
        
        % check for breaking waves
        if H/h > 0.78
            breaking(1,i) = 1;
        else
            breaking(1,i) = 0;
        end
        
        % check for weird waves
        pks = findpeaks(wave_p);
        trfs = findpeaks(-wave_p);
        if length(pks) == 1 && length(trfs) == 1
            weird_wave(1,i) = 0;
        else
            weird_wave(1,i) = 1;
        end
        % calculate q based on parameterized inputs
        [q_param_vec(1,i),s,Omegatc] = qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta,ws,param,Omegatc);
        q_mean_param(1,i) = nanmean(q_param_vec);
        
        % extract parameterized u and t and store in cell array
        u_param{1,i} = s.uw;
        t_param{1,i} = s.t;
        
    else % if p and u are not correlated, everything is NaN
    end
end

% create vector of indices
x = 1:length(q_vec);

% identify transport due to breaking waves
x_breaking = find(breaking == 1);
q_breaking = q_vec(breaking == 1);

% identify transport due to weird waves
x_weird = find(weird_wave == 1);
q_weird = q_vec(weird_wave == 1);

q_combined = [q_vec', q_param_vec'];
%% test plot
close all
figure()
plot(q_vec,'LineWidth',2)
hold on
plot(q_param_vec,'LineWidth',2)
hold off
yline(0,'--k')
title('Calculated transport per wave - Duck 94 Onshore Migration')
xlabel('wave')
ylabel('transport (m^2/s)')
set(gca,'FontSize',18)
legend('field data','parameterized')

figure()
bar(q_combined)
xlim([350 400])
title('Calculated transport per wave - Duck 94 Onshore Migration')
xlabel('wave')
ylabel('transport (m^2/s)')
set(gca,'FontSize',18)
legend('field data','parameterized')

figure()
histogram(abs(q_vec))
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Magnitude of transport (log log scale) - Duck 94 Onshore Migration')
xlabel('transport (m^2/s)')
ylabel('frequency')
set(gca,'FontSize',18)

figure()
plot(q_mean,'LineWidth',2)
hold on
plot(q_mean_param,'LineWidth',2)
hold off
title('cumulative mean transport - Duck 94 Onshore Migration')
ylabel('q mean(m^2/s)')
xlabel('wave')
yline(0,'--k')
set(gca,'FontSize',18)
legend('field data','parameterized')

%%
% plots a set of chosen waves
k = [941 624 91 90 660 31 100 905];
r = ceil(length(k)/2);

figure()
tiledlayout(r,2)
for j = 1:length(k)
    i = k(j);
    
    wave_t = t(upzc(i):upzc(i+1));
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1));

    nexttile
    plot(wave_t,wave_u,'k','LineWidth',2)
    hold on
    plot(wave_t,wave_p,'--k','LineWidth',2)
    hold off
    yline(0)
    title(num2str(i))
    if j==1
        legend('velocity', 'pressure')
    else
    end
end
