clear
close all

b0926 = matfile('bkgd0181_0926_1600EST.mat');
addpath('/Users/jack/Documents/GitHub/1DVarS/src/Uwave_ruessink2012')
addpath('/Users/jack/Documents/GitHub/1DVarS/src/qtransModels/vanderA2013')
addpath('/Users/jack/Documents/GitHub/1DVarS/src')

x_0926 = b0926.x;
d50_0926 = b0926.d50;
d90_0926 = b0926.d90;
h_0926 = b0926.h;
tanbeta_0926 = b0926.tanbeta;
Hrms_0926 = b0926.Hrms;
kabs_0926 = b0926.kabs;
omega_0926 = b0926.omega;
udelta_0926 = b0926.udelta_w;
delta_0926 = b0926.delta_bl;
ws_0926 = b0926.ws;
Aw_0926 = b0926.Aw;
Sw_0926 = b0926.Sw;
Uw_0926 = b0926.Uw;
Q_0926 = b0926.Q;
Q0_0926 = b0926.Q0;

figure()
tiledlayout(3,1)

nexttile
plot(x_0926,Q_0926)
title('Q')
xlabel('position (m)')
ylabel('Q (m^2/s')

nexttile
plot(x_0926, -h_0926)
title('h')

nexttile
plot(x_0926,Hrms_0926)
title('Hrms')
%%
clear distQ
clear qs
clear tl_qs
clear dHdQ

x1 = 225;
x2 = 350;
dx = 5;

xvec = x1:dx:x2;
idxmax = (x2-x1)/5 + 1;

for idx = 1:idxmax
    loc = x1 + (idx-1)*dx;
    locvec(1,idx) = loc; 


    pos_0926 = find(x_0926 == loc);
    Q0 = Q0_0926(pos_0926);

    param.alpha = 8.2; 
    param.xi = 1.7;
    param.m = 11;
    param.n = 1.2;
    param.nosusp = 1;

    d50 = d50_0926(pos_0926);
    d90 = d90_0926(pos_0926);
    h = h_0926(pos_0926);
    tanbeta = tanbeta_0926(pos_0926);
    Hrms = Hrms_0926(pos_0926);
    kabs = kabs_0926(pos_0926);
    k = kabs;
    omega = omega_0926;
    udelta = udelta_0926(pos_0926,:);
    delta = delta_0926(pos_0926);
    ws = ws_0926(pos_0926);
    B = Hrms/(sqrt(2));

    tl_d50 = 0;
    tl_d90 = 0;
    tl_h = 0;
    tl_tanbeta = 0;
    tl_Hrms = 1;
    tl_kabs = 0;
    tl_k = 0;
    tl_omega = 0;
    tl_udelta = [0,0];
    tl_delta = 0;
    tl_ws = 0;
    tl_Aw = 0;
    tl_Sw = 0;
    tl_Uw = 0;
    tl_param.alpha = 0;
    tl_param.xi = 0;
    tl_param.m = 0;
    tl_param.n = 0;
    tl_Hmo = 1;

    H = (0.1:0.1:2.5);

    for i = 1:length(H)
        Hmo = H(i);
        Hrms = H(i);
        [Aw,Sw,Uw,workspcU] = Uwave_ruessink2012_params(Hmo,k,omega,h);
        [qs(idx,i),workspc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);
        bkgdU = workspcU;
        [tl_Aw,tl_Sw,tl_Uw]=tl_Uwave_ruessink2012_params(tl_Hmo,tl_k,tl_omega,tl_h,bkgdU);
        bkgd = workspc;
        tl_qs(idx,i) = tl_qtrans_vanderA(tl_d50,tl_d90,tl_h,tl_tanbeta,tl_Hrms,tl_kabs,tl_omega,tl_udelta,tl_delta,tl_ws,tl_Aw,tl_Sw,tl_Uw,tl_param,bkgd);
        dHdQ(idx,i) = 1/tl_qs(idx,i);
        rayfcn = (Hmo/B^2)*exp((-Hmo^2)/(2*B^2));
        distQ(idx,i) = abs(dHdQ(idx,i)).*rayfcn;
%         distQ(distQ==inf) = [];
    end
    
    % linear interpolation for infinite tail
    infidx = find(distQ(idx,:) == Inf);
    infidxend = max(infidx);
    distQx1 = infidxend + 1;
    distQy1 = distQ(idx,distQx1);
    distQx2 = infidxend + 2;
    distQy2 = distQ(idx,distQx2);
    distQm = (distQy2-distQy1)/(distQx2-distQx1);
    
    for i = 1:length(infidx)
        distQlin(idx,infidx(i)) = distQm*(infidx(i) - distQx1) + distQy1;
    end
    for i = distQx1:length(H)
        distQlin(idx,i) = distQ(idx,i);
    end
    
    % exponential fit
    if ~isempty(infidx) 
        f = fit(H(distQx1+2:distQx1+7)',distQ(idx,distQx1+1:distQx1+6)','exp1');
        for i = 1:length(infidx)
            distQexp(idx,infidx(i)) = f.a*exp(f.b*H(infidx(i)));
        end
        for i = distQx1:length(H)
            distQexp(idx,i) = distQ(idx,i);
        end
    else
        distQexp(idx,:) = distQ(idx,:);
    end
    
    
end 

figure()
hold on
for idx = 1:idxmax
    plot(H,distQ(idx,:),'-k')
    plot(H,distQexp(idx,:),':b')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    integrate(1,idx) = trapz(qs(idx,:),distQexp(idx,:));
end
hold off


% sample plot x = 270m
figure()
tiledlayout(4,1)
nexttile
plot(H,qs(10,:))
set(gca,'xscale','log')
set(gca,'yscale','log')
title('transport x = 270m')
nexttile
plot(H,tl_qs(10,:))
set(gca,'xscale','log')
set(gca,'yscale','log')
title('dq/dH')
nexttile
plot(H,dHdQ(10,:))
title('dH/dq')
set(gca,'yscale','log')
set(gca,'xscale','log')
nexttile
plot(H,distQ(10,:))
xlabel('H (m)')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('distribution of Q')

% sample plot x = 225
figure()
tiledlayout(4,1)
nexttile
plot(H,qs(1,:))

title('transport x = 225m')
nexttile
plot(H,tl_qs(1,:))
set(gca,'xscale','log')
set(gca,'yscale','log')
title('dq/dH')
nexttile
plot(H,dHdQ(1,:))
title('dH/dq')
set(gca,'yscale','log')
set(gca,'xscale','log')
nexttile
plot(H,distQ(1,:))
xlabel('H (m)')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('distribution of Q')

% qbar = trapz(qs,qs.*distQ);
% qvar = trapz(qs,distQ.*(qs-qbar).^2);
% qstd = sqrt(qvar);
    
%     qbar(1,idx) = trapz(qs(idx,:),qs(idx,:).*distQ(idx,:));
%     qvar(1,idx) = trapz(qs(idx,:),distQ(idx,:).*(qs(idx,:)-qbar(1,idx)).^2);
%     qstd(1,idx) = sqrt(qvar(1,idx));
% end



% figure()
% plot(locvec,qbar)
% xlim([min(locvec) max(locvec)])
% title('qbar')
% xlabel('position (m)')
% 
% figure()
% plot(locvec,qvar)
% xlim([min(locvec) max(locvec)])
% title('qvar')
% xlabel('position (m)')

% figure()
% plot(qs,distQ)
% title('distribution of Q')
% set(gca,'xscale','log')
% set(gca,'yscale','log')