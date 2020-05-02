%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.4 - Identification of AR processes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 300]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);


%% -----------------------------------------------------------------------%
% 4.4.1 - Implement the given structure


%Synthesis
    N=1000;
    h=randn(N,1);

    % AR model coefficients
    a=[1 0.9 0.2];
    b=1;

    % Filter h WN 
    x=filter(b,a,h);

 % ANALYSIS
p=length(a)-1; %ignore 1
mu=0.01;

[y, error , w_adaptive] = adaptive_lms(x,mu,p);


Lx=length(x);
n=[0:Lx-1]; %discrete time axis

figure(1)
plot(n,w_adaptive,'Linewidth',1.5); hold on;
yline(-a(3),'k--')
yline(-a(2),'k--')
xlabel('Discrete time')
ylabel('Magnitude')
title(['Coefficient Evolution over time using constant \mu = ' num2str(mu)])
legend('\alpha_{1}','\alpha_{2}');
legend('Fontsize',13);
% savefig(figure(1),'figures/fig4_5.fig')
% saveas(figure(1),'figures/forlatex/fig4_5','epsc')


%% -----------------------------------------------------------------------%
% 4.4.22 - effect of mu

mu_tested=[0.01 0.04 0.08 0.12];
Ntests=length(mu_tested);
figure(2)
for m=1:Ntests
    subplot(2,2,m)
    mu=mu_tested(m);
    [y, error , w_adaptive] = adaptive_lms(x,mu,p);

    plot(n,w_adaptive,'Linewidth',1.5); hold on;
    yline(-a(3),'k--')
    yline(-a(2),'k--')
    xlabel('Discrete time')
    ylabel('Magnitude')
    title(['AR coefficient evolution using \mu = ' num2str(mu)])
    
end
legend('\alpha_{1}','\alpha_{2}','orientation','horizontal','Fontsize',12);
legend('box','off')

savefig(figure(2),'figures/fig4_6.fig')
% %% critical values of mu
% p=2;
%     [acfx, lags_xx]=xcorr(x,p,'unbiased');
% 
%     [acf_xz, lags_xz] = xcorr(x,h,p,'unbiased');
%     pzx=flip(acf_xz(lags_xz<=0));
% 
%     Rxx=zeros(p+1);
%     for column=0:-1:-p
%        Rxx(:,-column+1)=acfx(lags_xx>=column & lags_xx<=p+column);
%     end
% E = eig(Rxx);
% mean_convergence_crit_mu=2/max(E)
% ms_convergence_crit_mu=2/2*var(