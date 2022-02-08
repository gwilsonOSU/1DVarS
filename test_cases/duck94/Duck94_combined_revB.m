clear
close all

b0926 = matfile('bkgd0181_0926_1600EST.mat');
b1003 = matfile('bkgd0281_1003_2200EST.mat');

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

x_1003 = b1003.x;
d50_1003 = b1003.d50;
d90_1003 = b1003.d90;
h_1003 = b1003.h;
tanbeta_1003 = b1003.tanbeta;
Hrms_1003 = b1003.Hrms;
kabs_1003 = b1003.kabs;
omega_1003 = b1003.omega;
udelta_1003 = b1003.udelta_w;
delta_1003 = b1003.delta_bl;
ws_1003 = b1003.ws;
Aw_1003 = b1003.Aw;
Sw_1003 = b1003.Sw;
Uw_1003 = b1003.Uw;
Q_1003 = b1003.Q;
Q0_1003 = b1003.Q0;

%
close all

pos_0926 = find(x_0926 == 335);
pos_1003 = find(x_1003 == 335);
Q0 = [Q0_0926(pos_0926), Q0_1003(pos_1003)];

param.alpha = 8.2; 
param.xi = 1.7;
param.m = 11;
param.n = 1.2;
param.nosusp = 1;

N = 10;
n = 4000;

init = {nan(N,n),nan(N,n)};
q_re = init;
q_mean = init;
q_bulk = [nan,nan];

t_est = 2*N*n/2000*6/60;
prompt = ['Estimated runtime is ',num2str(t_est), ' minutes. Continue? [Y]'];
str = input(prompt,'s');

if str == 'Y'
    
    c = fix(clock);
    disp(['date: ',num2str(c(2)),'/',num2str(c(3))])
    disp(['time: ',num2str(c(4)),':',num2str(c(5))])
    tic
    
    for o = 1:2
        
        if o == 1

                d50 = d50_0926(pos_0926);
                d90 = d90_0926(pos_0926);
                h = h_0926(pos_0926) + b0926.tide;
                tanbeta = tanbeta_0926(pos_0926);
                kabs = kabs_0926(pos_0926);
                omega = omega_0926;
                udelta = udelta_0926(pos_0926,:);
                delta = delta_0926(pos_0926);
                ws = ws_0926(pos_0926);
                Hrms = Hrms_0926(pos_0926);

        elseif o == 2

                d50 = d50_1003(pos_1003);
                d90 = d90_1003(pos_1003);
                h = h_1003(pos_1003) + b1003.tide;
                tanbeta = tanbeta_1003(pos_1003);
                kabs = kabs_1003(pos_1003);
                omega = omega_1003;
                udelta = udelta_1003(pos_1003,:);
                delta = delta_1003(pos_1003);
                ws = ws_1003(pos_1003);
                Hrms = Hrms_1003(pos_1003);
                
        else

        end

        parfor i = 1:N

            Omegatc = 0;

            

            k = kabs;
            B = Hrms/(sqrt(2));
            Hmo = 1.4*Hrms;
            [Aw,Sw,Uw,~] = Uwave_ruessink2012_params(Hmo,k,omega,h);
            q_bulk(o) = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);
            
%             for ai = 1:5
                
%                 for bi = 1:5
                    
                    for j = 1:n

                        H = raylrnd(B);
                        
                        if H/h >= 0.78
                            
                            q_re{o}(i,j) = nan;
                            
                        else
                            
%                         a = 1.15;
%                         b = 0.85;

%                         a = 1.05 + ai./20;
%                         b = 0.75 + bi./20;

                        Hrms = H;
                        Hmo = H;

                        [Aw,Sw,Uw,~] = Uwave_ruessink2012_params(Hmo,k,omega,h);
%                         [q_re{o}(i,j),workspace_random] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param);
                        [q_re{o}(i,j),workspace_random,Omegatc] = qtrans_vanderA(d50,d90,h,tanbeta,Hrms,kabs,omega,udelta,delta,ws,Aw,Sw,Uw,param,Omegatc);
                        q_mean{o}(i,j) = nanmean(q_re{o}(i,1:j));
                        
                        end
                    end
                    
%                     q_abmean{o}(ai,bi) = mean(q_re{o}(1,:));

%                 end
                
%             end
           

        end
        
%         ab{o} = abs((q_abmean{o}-q_bulk(o))/q_bulk(o));
%         abidx(o) = find(ab{o} == min(min(ab{o})));
        
    end
    
%     ab_combined = ab{1} + ab{2};
%     ab_combinedidx = find(ab_combined == min(min(ab_combined)));
    
    toc
    
else
    
end



for o = 1:2
    
    figure()
    
    if o==1
        txt = 'onshore';
    else
        txt = 'offshore';
    end
    
    for i = 1:N
        
        plot(q_mean{o}(i,:),'LineWidth',2,'color',[0.5 0.5 0.5])
%         yline(q_abmean{o},'--r')
        hold on
    
    end
    
    hold off
    
    yline(Q0(o),'--k','LineWidth',2)
    yline(q_bulk(o),'--b')
    
    title('cumulative average transport')
    subtitle(txt)
    ylabel('q mean(m^2/s)')
    xlabel('wave')
%     ylim([-1E-4 1E-4])
    
end

figure()
histogram(q_re{1})
set(gca,'yscale','log')
set(gca,'xscale','log')
