clear 
close all

% paths
addpath(genpath('/home/shusin8/users/wilsongr/Duck94_SteveH/'))
addpath(genpath('/home/server/student/homes/aldricja/1DVarS'))

% load files
m_on = matfile('0926.mat');
dat_on = m_on.dat;
pos_on = m_on.pos;
m_off = matfile('1003.mat');
dat_off = m_off.dat;
pos_off = m_off.pos;

%%
fs = 2;

% extract data from files
time_on = 6; 
time_off = 4; 
loc_on = 10;
loc_off = 8;

% check u and p sensor positions for both cases
u_pos_on = pos_on(time_on,2,loc_on,1);
p_pos_on = pos_on(time_on,1,loc_on,1);
u_pos_off = pos_off(time_off,2,loc_off,1);
p_pos_off = pos_off(time_off,1,loc_off,1);

% time series for onshore
u_on = dat_on(time_on,2,loc_on,:);
u_on = u_on(:);
p_on = dat_on(time_on,1,loc_on,:);
p_on = p_on(:);
t_on = 0:1/fs:(length(u_on)-1)/2;
t_on = t_on';
bed_depth_on = pos_on(time_on,1,loc_on,4);
p_depth_on = pos_on(time_on,1,loc_on,3);
h_offset(1,1) = bed_depth_on - p_depth_on;

% time series for offshore
u_off = dat_off(time_off,2,loc_off,:);
u_off = u_off(:);
p_off = dat_off(time_off,1,loc_off,:);
p_off = p_off(:);
t_off = 0:1/fs:(length(u_off)-1)/2;
t_off = t_off';
bed_depth_off = pos_off(time_off,1,loc_off,4);
p_depth_off = pos_off(time_off,1,loc_off,3);
h_offset(1,2) = bed_depth_off - p_depth_off;

% combine data into cell arrays
u{1} = u_on;
u{2} = u_off;
p{1} = p_on;
pav(1,1) = mean(p{1});
p{2} = p_off;
pav(1,2) = mean(p{2});
t{1} = t_on;
t{2} = t_off;

% inputs that are the same for all cases
hav_vec = pav + h_offset;
d50 = 0.00018; 
d90 = 0.00024;
tanbeta_vec = [0.0276, 0.0227]; 
g = 9.81; 
rho = 1000; 
ws = ws_brownLawler(0.8*d50); 
wsc = ws;
wst = ws;
delta_on = .1262;
delta_off = 0.0662;
delta_vec = [delta_on,delta_off];
param.alpha = 8.2; 
param.xi = 1.7;
param.m = 11;
param.n = 1.2;
eta = nan;
lambda = nan;

cell_init = {nan,nan};
u_mean = cell_init;
p_mean = cell_init;
uWB = cell_init;
uLP = cell_init;
pWB = cell_init;
correlated = cell_init;
T_vec = cell_init;
H_vec = cell_init;
for o = 1:2
    tanbeta = tanbeta_vec(o);
    Omegatc = 0;
    % filter u and p
    % define moving mean for u and p
    filter_time = 512;
    p_mean{o} = movmean(p{o}, filter_time);

    % enter wave band filter parameters
    fcA = 1/20;  % wave band lower bound in Hz
    fcB = 0.5;   % wave band upper bound in Hz
    [bA,aA] = butter(5,fcA/(fs/2),'high');
    [bB,aB] = butter(5,fcB/(fs/2),'high');
    [bC,aC] = butter(5,fcB/(fs/2),'low');
    
    % detrend and filter u
    % wave band filter u 
    x = detrend(u{o});
    xA = filtfilt(bA,aA,x);
    uWB{o} = xA - filtfilt(bB,aB,xA);

    % low pass filter u
    x = detrend(u{o});
    trend = u{o} - detrend(u{o});
    xA = filtfilt(bC,aC,x);
    uLP{o} = xA + trend;

    % detrend and filter p
    x = detrend(p{o});
    xA = filtfilt(bA,aA,x);
    pWB{o} = xA - filtfilt(bB,aB,xA);
    
    % check correlation before calculating q
    % Chop into waves roughly based on p zero crossings
    upzc = find(diff(sign(pWB{o}))>0);
    init = NaN(1,length(upzc)-1);
    correlated{o} = init;
    % loop to check correlation
    
    for i = 1:(length(upzc)-1)
        
        % grabbing rough t, u, and p for each wave
        wave_t = t{o}(upzc(i):upzc(i+1));
        wave_uWB = uWB{o}(upzc(i):upzc(i+1));
        wave_u = uLP{o}(upzc(i):upzc(i+1));
        wave_p = pWB{o}(upzc(i):upzc(i+1)); 

        % checking correlation between p and u and calculating other inputs
        r = corrcoef(wave_uWB,wave_p);
        r = abs(r(2,1));
        
        if r >= 0.9
            
            % calculate T and store as cell array
            x1 = t{o}(upzc(i));
            x2 = t{o}(upzc(i)+1);
            y1 = pWB{o}(upzc(i));
            y2 = pWB{o}(upzc(i)+1);
            m = (y2-y1)/(x2-x1);
            t00 = x1 - y1/m; 
            
            x1 = t{o}(upzc(i+1));
            x2 = t{o}(upzc(i+1)+1);
            y1 = pWB{o}(upzc(i+1));
            y2 = pWB{o}(upzc(i+1)+1);
            m = (y2-y1)/(x2-x1);
            t02 = x1 - y1/m;
            
            T = t02 - t00;
            T_vec{o}(1,i) = T;
            
            % calculate H and store as cell array
            H = max(wave_p) - min(wave_p);
            H_vec{o}(1,i) = H;
            
            udelta = [mean(wave_u(2:end)),0];
            uw = wave_u - udelta(1);
            delta = delta_vec(o);
            
            % uhat, uhatc, and uhatt based on definitions
            uhat = sqrt(2*mean(uw.^2));        
            uhatc = max(uw);        
            uhatt = abs(min(uw));
            
            % interpolating t for mid point of wave
            midzc = find(diff(sign(wave_p))<0);
            x1 = wave_t(midzc);
            x2 = wave_t(midzc+1);
            y1 = wave_p(midzc);
            y2 = wave_p(midzc+1);
            m = (y2-y1)/(x2-x1);
            t01 = x1 - y1/m;

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
            h = p_mean{o}(upzc(i)) + h_offset(o);
            omega = 2*pi/T;
            k = fzero( @(k) omega^2 - g*k*tanh(k*h), 0.1);
            kabs = k;
            c = omega/k;
            
            % calculate q and store as cell array
            [q_vec{o}(1,i),~,Omegatc] = qtrans_vanderA_onewave(uw,d50,d90,wsc,wst,udelta,delta,uhat,uhatc,uhatt,T,Tc,Tt,Ttu,Tcu,c,eta,lambda,tanbeta,param,Omegatc);
            q_mean{o}(1,i) = nanmean(q_vec{o});
            
            
        else
            T_vec{o}(1,i) = nan;   
            H_vec{o}(1,i) = nan;
            q_vec{o}(1,i) = nan;
            q_mean{o}(1,i) = nan;
        end
        
    end
end
%%
Tav_vec = [nanmean(T_vec{1}),nanmean(T_vec{2})];
q_Re = cell_init;
q_Re_mean = cell_init;
H_ray = cell_init;
N = 10;
n = 3.7E3;
time_est = N*n*2*12*218/533/1E5;
prompt = ['Estimated runtime is ',num2str(time_est), ' minutes. Continue? [Y]'];
str = input(prompt,'s');
if str == 'Y'
    tic
    % loop over offshore and onshore cases
    for o = 1:2
        tanbeta = tanbeta_vec(o);
        delta = delta_vec(o);
        omega = 2*pi/Tav_vec(o);
        h = hav_vec(o); 
        kabs = fzero( @(kabs) omega^2 - g*kabs*tanh(kabs*h), 0.1); 
        k = kabs;
        c = omega/kabs; 

        H_vec{o}(isnan(H_vec{o}))=[];
        ray = fitdist(H_vec{o}.','Rayleigh');
        
        Omegatc = 0;
        Hrms = sqrt(mean(H_vec{o}.^2,2));
        Hmo = 1.4*Hrms;
        E = 0.125*rho*g*Hrms^2; 
        undertow = -E/(rho*c*h);
        udelta = [undertow,0];
        [Aw,Sw,Uw,~]=Uwave_ruessink2012_params(Hmo,k,omega,h);
        [q_bulk(1,o),~,~] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);
        % inputs/outputs unique to each random wave
        for i = 1:N
            Omegatc = 0; 
            for j = 1:n
                Hrms = random(ray); 
                Hmo = 1.4*Hrms;
                H_ray{o}(i,j) = Hrms;
                E = 0.125*rho*g*Hrms^2; 
                undertow = -E/(rho*c*h);
                udelta = [undertow,0]; 
                [Aw,Sw,Uw,~]=Uwave_ruessink2012_params(Hmo,k,omega,h);
                [q_Re{o}(i,j),~,Omegatc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);
                q_Re_mean{o}(i,j) = nanmean(q_Re{o}(i,1:j));
            end
        end

        figure()
        if o==1
            txt = 'onshore';
        else
            txt = 'offshore';
        end
        plot(q_mean{o})
        hold on
        for i = 1:N
            plot(q_Re_mean{o}(i,:),'LineWidth',2,'color',[.5 .5 .5])
            hold on
        end
        hold off
        title('running average transport')
        subtitle(txt)
        ylabel('q mean(m^2/s)')
        xlabel('wave')
        yline(0,'--k')
        set(gca,'FontSize',18)


    end
    toc
else
end


