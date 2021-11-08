clear
close all

addpath(genpath('/Users/jack/Documents/Github/1DVarS'))

% Load p,t, uE, uN
uE = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','currentEast');
uN = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','currentNorth');
t = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','time');
p = ncread('FRF-ocean_WaveShape_Vector_Crab20161118_01.nc','pressure');

% Normalize t
t = t-t(1);

% Calculate u
u = -uE;
% theta = 71.8;
% u = uE*cosd(90-theta) + uN*sind(theta);
% v = uN*cosd(90-theta) - uE*sind(theta);

% define moving mean for u and p
filter_time = 512;
u_mean = movmean(u,filter_time);
p_mean = movmean(p, filter_time);

% enter wave band filter parameters
fs = 2;
fcA = 1/15;  % wave band lower bound in Hz
fcB = 0.7;   % wave band upper bound in Hz
[bA,aA] = butter(5,fcA/(fs/2),'high');
[bB,aB] = butter(5,fcB/(fs/2),'high');

% detrend and filter u
x = detrend(u);
xA = filtfilt(bA,aA,x);
uWB = xA - filtfilt(bB,aB,xA);

% detrend and filter p
x = detrend(p);
xA = filtfilt(bA,aA,x);
pWB = xA - filtfilt(bB,aB,xA);

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
init = zeros(1,length(upzc)-1);
uhat_vec = init;
uhatc_vec = init;
uhatt_vec = init;
T_vec = init;
Tc_vec = init;
Tt_vec = init;
Ttu_vec = init;
Tcu_vec = init;
q_vec = init;
H_vec = init;
breaking = init;
q_param_vec = init;
uhat_param_vec = init;
uhatc_param_vec = init;
uhatt_param_vec = init;
Tc_param_vec = init;
Tt_param_vec = init;
Tcu_param_vec = init;
Ttu_param_vec = init;
h_vec = init;
k_vec = init;
omega_vec = init;
c_vec = init;

cell_init = cell(1,length(upzc)-1);
u_param = cell_init;
t_param = cell_init;
udelta_cell = cell_init;

% loop to calculate other inputs and q
for i = 1:(length(upzc)-1)
    % grabbing rough t, u, and p for each wave
    wave_t = t(upzc(i):upzc(i+1));
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1)); 
    
    % checking correlation between p and u and calculating other inputs
    r = corrcoef(wave_u,wave_p);
    r = abs(r(2,1));
    if r >= 0.9
        % udelta as mean of u at rough start of wave 
        udelta = [u_mean(upzc(i)),0];
        udelta_cell{1,i} = udelta;
        
        % uhat, uhatc, and uhatt based on definitions
        uhat = sqrt(2*mean(wave_u.^2));
        uhat_vec(1,i) = uhat;
        
        uhatc = max(wave_u);
        uhatc_vec(1,i) = uhatc;
        
        uhatt = abs(min(wave_u));
        uhatt_vec(1,i) = uhatt;
        
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
        Tc_vec(1,i) = Tc;
        
        Tt = t02 - t01;
        Tt_vec(1,i) = Tt;
        
        % calculating Ttu and Tcu, with precautions to keep Ttu positive
        zcross = floor(t01);
        trough = find(wave_t>=zcross);
        [maxcrest,maxcrestidx] = max(wave_p);
        [mintrough,mintroughidx] = min(wave_p(trough));
        
        Ttu = wave_t(mintroughidx+length(wave_t)-length(trough)) - t01;
        Ttu_vec(1,i) = Ttu;
        
        Tcu = wave_t(maxcrestidx) - t00;
        Tcu_vec(1,i) = Tcu;
        
        % dispersion solver for c
        h = p_mean(upzc(i));
        h_vec(1,i) = h;
        
        omega = 2*pi/T;
        omega_vec(1,i) = omega;
        k = fzero( @(k) omega^2 - g*k*tanh(k*h), 0.1);
        k_vec(1,i) = k;
        kabs = k;
        c = omega/k;
        c_vec(1,i) = c;
        
        % calculate q
        q_vec(1,i) = qtrans_vanderA_onewave_ja(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param);
        
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
        
        % calculate q based on parameterized inputs
        [q_param_vec(1,i),s] = qtrans_vanderA(d50,d90,h,Hrms,kabs,omega,udelta,ws,param);
        
        % extract parameterized inputs and store in arrary
        uhat_param_vec(1,i) = s.uhat;
        uhatc_param_vec(1,i) = s.uhatc;
        uhatt_param_vec(1,i) = s.uhatt;
        Tc_param_vec(1,i) = s.Tc;
        Tt_param_vec(1,i) = s.Tt;
        Tcu_param_vec(1,i) = s.Tcu;
        Ttu_param_vec(1,i) = s.Ttu;
        
        % extract parameterized u and t and store in cell array
        u_param{1,i} = s.uw;
        t_param{1,i} = s.t;
        
    else % if p and u are not correlated, everything is NaN
        uhat_vec(1,i) = NaN;
        uhatc_vec(1,i) = NaN;
        uhatt_vec(1,i) = NaN;
        T_vec(1,i) = NaN;
        Tc_vec(1,i) = NaN;
        Tt_vec(1,i) = NaN;
        Ttu_vec(1,i) = NaN;
        Tcu_vec(1,i) = NaN;
        q_vec(1,i) = NaN;
        H_vec(1,i) = NaN;
        breaking(1,i) = NaN;
        q_param_vec(1,i) = NaN;
        uhat_param_vec(1,i) = NaN;
        uhatc_param_vec(1,i) = NaN;
        uhatt_param_vec(1,i) = NaN;
        Tc_param_vec(1,i) = NaN;
        Tt_param_vec(1,i) = NaN;
        Tcu_param_vec(1,i) = NaN;
        Ttu_param_vec(1,i) = NaN;
        u_param{1,i} = NaN;
        t_param{1,i} = NaN;
        h_vec(1,i) = NaN;
        k_vec(1,i) = NaN;
        omega_vec(1,i) = NaN;
        c_vec(1,i) = NaN;
        udelta_cell{1,i} = [NaN NaN];
    end
end

% create vector of indices
x = 1:length(q_vec);

% identify transport due to breaking waves
x_breaking = find(breaking == 1);
q_breaking = q_vec(breaking == 1);

% plot transport
figure(1)
plot(q_vec,'-ok')
hold on
plot(q_param_vec,'-ob')

% flag breaking
scatter(x_breaking,q_breaking,'r','filled')
hold off

% plot formatting
title('Calculated transport per wave')
xlabel('wave identifier')
ylabel('transport (m^2/s)')
yline(0)
legend('inputs from data','parameterized inputs','breaking waves')

%%
% scatter q vs q
figure(2)
scatter(q_vec,q_param_vec,'k','filled')
limit = max(abs([q_vec q_param_vec]));
ylim([-limit limit])
xlim([-limit limit])
yline(0,'--')
xline(0,'--')
title('transport calculated from parameterization vs transport calculated from field data')
xlabel('transport based on inputs from field data (m^2/s)')
ylabel('transport based on parameterized inputs (m^2/s)')
x_onetoone = linspace(-limit,limit);
y_onetoone = x_onetoone;
hold on
plot(x_onetoone,y_onetoone,'--k')
hold off

%%
% scatter data inputs v parameterized inputs
data_inputs = [Tc_vec; uhat_vec; Tt_vec; uhatc_vec; Tcu_vec; uhatt_vec; Ttu_vec];
param_inputs = [Tc_param_vec; uhat_param_vec; Tt_param_vec; uhatc_param_vec; Tcu_param_vec; uhatt_param_vec; Ttu_param_vec;];
title_array = ["Tc" "uhat" "Tt" "uhatc" "Tcu" "uhatt" "Ttu"];

figure(3)
tiledlayout(4,2)

for i = 1:7
    nexttile
    scatter(data_inputs(i,:),param_inputs(i,:),'k','filled')
    title(title_array(i))
    limit = max(abs([data_inputs(i,:) param_inputs(i,:)]));
    ylim([0 limit])
    xlim([0 limit])
    yline(0,'--')
    xline(0,'--')
    xlabel('inputs from data')
    ylabel('parameterized inputs')
    x_onetoone = linspace(-limit,limit);
    y_onetoone = x_onetoone;
    hold on
    plot(x_onetoone,y_onetoone,'--k')
    hold off
end

%%
% plots a set of chosen waves
k = [14    29    64    92   134   162   202   222   230];
r = ceil(length(k)/2);

figure()
tiledlayout(r,2)
for j = 1:length(k)
    i = k(j);
    
    wave_t = t(upzc(i):upzc(i+1));
    wave_u = uWB(upzc(i):upzc(i+1));
    wave_p = pWB(upzc(i):upzc(i+1));

    % interpolating t for beginning of wave
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

    % interpolating t for mid point of wave
    midzc = find(diff(sign(wave_p))<0);
    x11 = wave_t(midzc);
    x21 = wave_t(midzc+1);
    y11 = wave_p(midzc);
    y21 = wave_p(midzc+1);
    m1 = (y21-y11)/(x21-x11);
    t01 = x11 - y11/m1;


    % calculating Ttu and Tcu, with precautions to keep Ttu positive
    zcross = floor(t01);
    trough = find(wave_t>=zcross);
    [maxcrest,maxcrestidx] = max(wave_p);
    [mintrough,mintroughidx] = min(wave_p(trough));

    nexttile
    plot(wave_t,wave_u,'k')
    hold on
    plot(wave_t,wave_p,'--k')
    plot(t_param{i}+wave_t(1),u_param{i},'b')
    ylim([-1 1.25])
    xlim([wave_t(1) wave_t(1)+12])
    hold off
    yline(0)
    xline(t01,'--')
    % xline(wave_t(1)+Tc_param_vec(i),'--b')
    xline(wave_t(maxcrestidx),':')
    % xline(wave_t(1)+Tcu_param_vec(i),':b')
    xline(wave_t(mintroughidx+length(wave_t)-length(trough)),'-.')
    % xline(wave_t(1)+Tc_param_vec(i)+Ttu_param_vec(i),'-.b')
    title(['wave ' num2str(i)])
    if j==1
        legend('real velocity', 'real pressure','parameterized velocity')
    else
    end
end

%%
% swaps inputs one at a time

% rearranging inputs for grouping
data_inputs_swap(1,:) = data_inputs(1,:);
data_inputs_swap(2,:) = data_inputs(3,:);
data_inputs_swap(3,:) = data_inputs(2,:);
data_inputs_swap(4,:) = data_inputs(4,:);
data_inputs_swap(5,:) = data_inputs(6,:);
data_inputs_swap(6,:) = data_inputs(5,:);
data_inputs_swap(7,:) = data_inputs(7,:);

param_inputs_swap(1,:) = param_inputs(1,:);
param_inputs_swap(2,:) = param_inputs(3,:);
param_inputs_swap(3,:) = param_inputs(2,:);
param_inputs_swap(4,:) = param_inputs(4,:);
param_inputs_swap(5,:) = param_inputs(6,:);
param_inputs_swap(6,:) = param_inputs(5,:);
param_inputs_swap(7,:) = param_inputs(7,:);

input_matrix = data_inputs_swap; % initializing w inputs from data

q_vec_inputswap = zeros(5,236);
for j = 1:5 
    if j <= 2
        input_matrix(j,:) = param_inputs_swap(j,:); % swaps Tc or Tt
        input_matrix(j+5,:) = param_inputs_swap(j+5,:); % swaps Tcu or Ttu
    else
        input_matrix(j,:) = param_inputs_swap(j,:);
    end
    for i = 1:length(upzc)-2 
        if ~isnan(c_vec(i)) == 1 % checking for NaN
            % these two inputs don't change
            udelta = udelta_cell{i};
            c = c_vec(i);
            
            % inputs that might be swapped:
            Tc = input_matrix(1,i);
            Tt = input_matrix(2,i);
            uhat = input_matrix(3,i);
            uhatc = input_matrix(4,i);
            uhatt = input_matrix(5,i);
            Tcu = input_matrix(6,i);
            Ttu = input_matrix(7,i);
            
            % calculate new q
            q_vec_inputswap(j,i) = qtrans_vanderA_onewave_ja(d50,d90,wsc,wst,udelta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,param);
        else
            q_vec_inputswap(j,i) = NaN;
        end
    end
    input_matrix = data_inputs_swap; % resets to data inputs after calculating q
end

title_array_swap = ["Tc and Tcu" "Tt and Ttu" "uhat" "uhatc" "uhatt"];

figure(5)
tiledlayout(3,2)

for i = 1:5
    nexttile
    scatter(q_vec(1:end-1),q_vec_inputswap(i,:),'k','filled')
    % limit = max(abs([q_vec(1:end-1) q_vec_inputswap(i,:)]));
    limit = 1E-4;
    ylim([-limit limit])
    xlim([-limit limit])
    yline(0,'--')
    xline(0,'--')
    title(title_array_swap(i))
    xlabel('q - all inputs from data')
    ylabel('q - one parameterized input')
    x_onetoone = linspace(-limit,limit);
    y_onetoone = x_onetoone;
    hold on
    plot(x_onetoone,y_onetoone,'--k')
    hold off
end

%%

figure(6)
plot(x(1:end-1),q_vec(1:end-1),'-o',x(1:end-1),q_vec_inputswap(1,:),'-o',x(1:end-1),q_vec_inputswap(2,:),'-o',...
    x(1:end-1),q_vec_inputswap(3,:),'-o',x(1:end-1),q_vec_inputswap(4,:),'-o',x(1:end-1),q_vec_inputswap(5,:),'-o')
legend(["data",title_array_swap])
title('net transport swapping one parameter at a time')
xlabel('wave identifier')
ylabel('transport (m^2/s)')
yline(0)