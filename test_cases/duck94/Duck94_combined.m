%% file loading stuff
clear 
close all

addpath(genpath('/home/shusin8/users/wilsongr/Duck94_SteveH/'))
addpath(genpath('/home/server/student/homes/aldricja/1DVarS'))

m_on = matfile('0926.mat');
dat_on = m_on.dat;
pos_on = m_on.pos;

m_off = matfile('1003.mat');
dat_off = m_off.dat;
pos_off = m_off.pos;

%% define variables
fs = 2;
time_on = 6; % time in 8hr blocks 
loc_on = 10;
u_pos = pos_on(time_on,2,loc_on,1);
p_pos = pos_on(time_on,1,loc_on,1);
u_on = dat_on(time_on,2,loc_on,:);
u_on = u_on(:);
p_on = dat_on(time_on,1,loc_on,:);
p_on = p_on(:);
t_on = 0:1/fs:(length(u_on)-1)/2;
t_on = t_on';

time_off = 4; % time in 8hr blocks 
loc_off = 8;
u_pos_off = pos_off(time_off,2,loc_off,1);
p_pos_off = pos_off(time_off,1,loc_off,1);
u_off = dat_off(time_off,2,loc_off,:);
u_off = u_off(:);
p_off = dat_off(time_off,1,loc_off,:);
p_off = p_off(:);
t_off = 0:1/fs:(length(u_off)-1)/2;
t_off = t_off';

%% combine data into cell arrays
u{1} = u_on;
u{2} = u_off;
p{1} = p_on;
p{2} = p_off;
t{1} = t_on;
t{2} = t_off;

%% check plot before filtering
close all
for i = 1:2
    if i==1
        txt = 'onshore';
    else
        txt = 'offshore';
    end
    p_plot = detrend(p{i});
    figure()
    plot(t{i},u{i},t{i},p_plot)
    legend('velocity','pressure')
    xlim([0 100])
    title('sanity check')
    subtitle(txt)
    yline(0,'--k')
end

%% filter u and p
cell_init = {nan,nan};
u_mean = cell_init;
p_mean = cell_init;
uWB = cell_init;
uLP = cell_init;
pWB = cell_init;

for i = 1:2
    % define moving mean for u and p
    filter_time = 512;
    u_mean{i} = movmean(u{i},filter_time);
    p_mean{i} = movmean(p{i}, filter_time);

    % enter wave band filter parameters
    fcA = 1/20;  % wave band lower bound in Hz
    fcB = 0.5;   % wave band upper bound in Hz
    [bA,aA] = butter(5,fcA/(fs/2),'high');
    [bB,aB] = butter(5,fcB/(fs/2),'high');
    [bC,aC] = butter(5,fcB/(fs/2),'low');
    % detrend and filter u
    % wave band filter u 
    x = detrend(u{i});
    xA = filtfilt(bA,aA,x);
    uWB{i} = xA - filtfilt(bB,aB,xA);

    % low pass filter u
    x = detrend(u{i});
    trend = u{i} - detrend(u{i});
    xA = filtfilt(bC,aC,x);
    uLP{i} = xA + trend;

    % detrend and filter p
    x = detrend(p{i});
    xA = filtfilt(bA,aA,x);
    pWB{i} = xA - filtfilt(bB,aB,xA);
end

%% check correlation before calculating q
correlated = cell_init;
for i = 1:2
    % Chop into waves roughly based on p zero crossings
    upzc = find(diff(sign(pWB{i}))>0);
    init = NaN(1,length(upzc)-1);
    correlated{i} = init;
    % loop to check correlation
    for ii = 1:(length(upzc)-1)
        % grabbing rough t, u, and p for each wave
        wave_t = t{i}(upzc(ii):upzc(ii+1));
        wave_u = uWB{i}(upzc(ii):upzc(ii+1));
        wave_p = pWB{i}(upzc(ii):upzc(ii+1)); 

        % checking correlation between p and u and calculating other inputs
        r = corrcoef(wave_u,wave_p);
        r = abs(r(2,1));
        if r >= 0.9
            correlated{i}(1,ii) = 1;
        else
            correlated{i}(1,ii) = 0;
        end
    end
end

%% calculate q using pressure whenever possible
% initializing vectors for other inputs
q_vec = cell_init;
H_vec = cell_init;
T_vec = cell_init;
breaking = cell_init;
q_mean = cell_init;
weird_wave = cell_init;
q_param_vec = cell_init;
q_mean_param = cell_init;
udelta_vec = cell_init;

for i = 1:2
    % Chop into waves roughly based on p zero crossings
    upzc = find(diff(sign(pWB{i}))>0);
    
    % initialize vectors within cells to be NaN
    init = NaN(1,length(upzc)-1);
    q_vec{i} = init;
    H_vec{i} = init;
    T_vec{i} = init;
    breaking{i} = init;
    q_mean{i} = init;
    weird_wave{i} = init;
    q_param_vec{i} = init;
    q_mean_param{i} = init;
    udelta_vec{i} = init;
    
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
    eta = NaN;
    lambda = NaN;
    if i == 1
        delta = 0.05*5;
    else
        delta = 0.05;
    end
    tanbeta = 0.05; % constate for now, also from Greg
    Omegatc = 0;
    % loop to calculate other inputs and q
    for ii = 1:(length(upzc)-1)
        % grabbing rough t, u, and p for each wave
        wave_t = t{i}(upzc(ii):upzc(ii+1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% possible bug!! should be t(upzc(i)+1) ???
        wave_u = uLP{i}(upzc(ii):upzc(ii+1));
        wave_p = pWB{i}(upzc(ii):upzc(ii+1)); 

        % checking correlation between p and u and calculating other inputs
        if correlated{i}(1,ii) == 1
            % udelta as mean of u at rough start of wave 
            udelta = [mean(wave_u(2:end)),0];
            udelta_vec{i}(1,ii) = mean(wave_u(2:end));
            uw = wave_u - udelta(1);

            % uhat, uhatc, and uhatt based on definitions
            uhat = sqrt(2*mean(uw.^2));        
            uhatc = max(uw);        
            uhatt = abs(min(uw));

            % interpolating t for beginning of wave
            x10 = t{i}(upzc(ii));
            x20 = t{i}(upzc(ii)+1);
            y10 = pWB{i}(upzc(ii));
            y20 = pWB{i}(upzc(ii)+1);
            m0 = (y20-y10)/(x20-x10);
            t00 = x10 - y10/m0; 

            % interpolating t for end of wave
            x12 = t{i}(upzc(ii+1));
            x22 = t{i}(upzc(ii+1)+1);
            y12 = pWB{i}(upzc(ii+1));
            y22 = pWB{i}(upzc(ii+1)+1);
            m2 = (y22-y12)/(x22-x12);
            t02 = x12 - y12/m2;

            % calculating T
            T = t02 - t00;
            T_vec{i}(1,ii) = T;

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
            h = p_mean{i}(upzc(ii));        
            omega = 2*pi/T;
            k = fzero( @(k) omega^2 - g*k*tanh(k*h), 0.1);
            kabs = k;
            c = omega/k;

            % calculate q
            [q_vec{i}(1,ii),~,Omegatc] = qtrans_vanderA_onewave(uw,d50,d90,wsc,wst,udelta,delta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,eta,lambda,tanbeta,param,Omegatc);
            q_mean{i}(1,ii) = nanmean(q_vec{i});
            % calculate wave height
            H = max(wave_p) - min(wave_p);
            H_vec{i}(1,ii) = H;
            Hrms = H;
            Hmo = H;
            % check for breaking waves
            if H/h > 0.78
                breaking{i}(1,ii) = 1;
            else
                breaking{i}(1,ii) = 0;
            end

            % check for weird waves
            pks = findpeaks(wave_p);
            trfs = findpeaks(-wave_p);
            if length(pks) == 1 && length(trfs) == 1
                weird_wave{i}(1,ii) = 0;
            else
                weird_wave{i}(1,ii) = 1;
            end
            % calculate q based on parameterized inputs
            [Aw,Sw,Uw,~]=Uwave_ruessink2012_params(Hmo,k,omega,h);
            [q_param_vec{i}(1,ii),s,Omegatc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);
            q_mean_param{i}(1,ii) = nanmean(q_param_vec{i});
            
            % extract parameterized u and t and store in cell array
            u_param{i,ii} = s.uw;
            t_param{i,ii} = s.t;

        else % if p and u are not correlated, everything is NaN
        end
    end

    % create vector of indices
    x = 1:length(q_vec);

    % identify transport due to breaking waves
    x_breaking = find(breaking{i} == 1);
    q_breaking = q_vec{i}(breaking{i} == 1);

    % identify transport due to weird waves
    x_weird = find(weird_wave{i} == 1);
    q_weird = q_vec{i}(weird_wave{i} == 1);

%     q_combined = [q_vec{i}', q_param_vec{i}'];

end

%% test plot
close all
for i = 1:2
    if i==1
        txt = 'onshore';
    else
        txt = 'offshore';
    end
    figure()
    plot(q_vec{i},'LineWidth',2)
    hold on
    plot(q_param_vec{i},'LineWidth',2)
    hold off
    yline(0,'--k')
    title('Calculated transport per wave - Duck 94')
    xlabel('wave')
    ylabel('transport (m^2/s)')
    set(gca,'FontSize',18)
    legend('field data','parameterized')
    subtitle(txt)
    
    figure()
    histogram(abs(q_vec{i}))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title('Magnitude of transport (log log scale) - Duck 94')
    xlabel('transport (m^2/s)')
    ylabel('frequency')
    set(gca,'FontSize',18)
    subtitle(txt)

    figure()
    plot(q_mean{i},'LineWidth',2)
    hold on
    plot(q_mean_param{i},'LineWidth',2)
    hold off
    title('cumulative mean transport - Duck 94')
    ylabel('q mean(m^2/s)')
    xlabel('wave')
    yline(0,'--k')
    set(gca,'FontSize',18)
    legend('field data','parameterized')
    subtitle(txt)

end

%% Rayleigh distributed random waves
% create pdf using built in function, removing NaNs

n = 1; % number of runs in ensemble
N = 1E3; % number of waves in each run: 1E5 = ~1 week, 12 minutes run time

% initialize q and q_mean
init = NaN(n,N); 
cell_init = {init,init};
q = cell_init;
q_ray_mean = cell_init;
breaking = cell_init;
for i = 1:2
    H_vec{i}(isnan(H_vec{i}))=[];
    ray = fitdist(H_vec{i}.','Rayleigh');

    % Other inputs used to calculate q
    T_av = nanmean(T_vec{i});
    omega = 2*pi/T_av; % definition
    kabs = fzero( @(kabs) omega^2 - g*kabs*tanh(kabs*h), 0.1); % dispersion solver
    udelta = [mean(u_on),0];
    Omegatc = 0; % initialized as zero for the first wave
    for j = 1:n
        for ii=1:N
            Hrms = random(ray); % choose Rayleigh distributed random wave height
            Hmo = Hrms;
            [Aw,Sw,Uw,~]=Uwave_ruessink2012_params(Hmo,k,omega,h);           
            [q{i}(1,ii),s,Omegatc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);
            q_ray_mean{i}(1,ii) = nanmean(q{i});
            if Hrms/h > 0.78 % check for depth limited breaking
                breaking{i}(1,ii) = 1;
            else
                breaking{i}(1,ii) = 0;
            end
        end
    end
end

%% plots 
close all

for j = 1:2
    figure()
    for i = 1:n
        plot(q_mean{j}(i,:),'LineWidth',2)
        hold on
    end
    hold off
    title('running average transport')
    ylabel('q mean(m^2/s)')
    xlabel('wave')
    yline(0,'--k')
    yline(q_mean{j}(n,N))
    yline(1.1*q_mean{j}(n,N),'--k')
    yline(0.9*q_mean{j}(n,N),'--k')
    set(gca,'FontSize',18)

    figure()
    tiledlayout(3,1)
    nexttile
    histogram(abs(q{j}),'FaceColor','r')
    title('q, linear scale')
    nexttile
    histogram(abs(q{j}),'FaceColor','r')
    set(gca,'yscale','log')
    title('q, log scale y axis')
    nexttile
    histogram(abs(q{j}),'FaceColor','r')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title('q, log log scale')
end