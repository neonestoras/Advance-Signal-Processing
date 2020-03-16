%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 Cramer-Rao Lower Bound   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 310]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);



%%  2.4.1 - AR model of NASDAQ financial %% 
%-------------------------------------------------------------------------%
load 'NASDAQ.mat';
closing_prices = NASDAQ.Close;
N=length(closing_prices);
time = NASDAQ.Date;

figure (1)
plot(time, closing_prices);

figure(2)
subplot(1,2,1)
[a_or,e_or,rc_or] = aryule(closing_prices,10);
pacf_cp_org=-rc_or; %partial acf equal to rc*-1

%standardise the series (zero mean and unit variance)
stand_cp = (closing_prices - mean(closing_prices))./std(closing_prices);
[a_stand,e_stand,rc_stand] = aryule(stand_cp,10); % a_stand is 
pacf_cp_stand=-rc_stand;

stem(pacf_cp_org,'b','Linewidth',1.4); hold on;
stem(pacf_cp_stand,'r', 'MarkerFaceColor','red','MarkerEdgeColor','none','Linewidth',1.3);
yline(1.96/sqrt(288), '--g','Linewidth',1.25);  % 95% confidence interval 
yline(-1.96/sqrt(288), '--g','Linewidth',1.25);

title('PACF of the original and standarised data model'); 
ylabel('Correlation');
xlabel('Time-Lag (k)');
legend('\fontsize{10}PACF original series','\fontsize{10}PACF standarised series','\fontsize{10}95% confidence interval')
legend('Location','northeast');
hold off;


subplot(1,2,2)

mdl=[];
aic=[];
aic_c=[];

for p= 1:10 %i=1 for a0, p: parameters iteration
    [A,e(p)] = aryule(stand_cp,p);
    mdl (:,p) = log10(e(p)) + ((p*log10(N))/N); % minimum description length criterion
    aic (:,p) = log10(e(p)) + (2*p)/N; %Akaike information criterion
    aic_c (:,p) =  aic(:,p) + (((2*p)*(p+1))/(N-p-1)); %Corrected AIC 
end


plot(mdl,'b','LineWidth' ,2.2); hold on;
plot(aic,'r','LineWidth' ,1.35);
plot(aic_c,'LineWidth' ,1.5);
plot(log10(e),'k--','LineWidth' ,1.5)
legend('\fontsize{10}MDL', '\fontsize{10}AIC','\fontsize{10}AIC corrected','\fontsize{10}Loss function','Location','northwest')
title('MDL, AIC and AICc for standarised data')
xlabel('Model order');xlim([1 10])
ylabel('Prediction Error'); ylim([-1.8 -0.8])

% savefig(figure(2),'figures/fig2_14.fig')
% saveas(figure(2),'figures/forlatex/fig2_14','epsc')

%%  2.4.c.i - Plotting the CLRB for two given parameters %% 
%-------------------------------------------------------------------------%

%original data
Rx= xcorr(closing_prices,'unbiased');
max_or=max(Rx);

%standarised data
Rx_stand=xcorr(stand_cp,'unbiased');
max_stand=max(Rx_stand);

N=1001:-50:1; %Number of data points(order is reversed in order to be desplayed in correct order on heatmap)
var_noise = 1:50:1001; %True variance of input driving noise

for N_i=1:length(N)
    
    for var_i = 1:length(var_noise)
        CRLB_var(N_i,var_i)=(2*var_noise(var_i)^2)/N(N_i);
        CRLB_a(N_i,var_i)=(var_noise(var_i))/(N(N_i)*max_stand); 
    end
end


%creating the heatmap
figure(3)
subplot(1,2,1)
heatmap(var_noise,N,log10(CRLB_var),'Colormap', parula(256),'ColorbarVisible','on','YLabel','N','XLabel','Noise Variance');
title('CRLB for estimate of variance of driving noise')
subplot(1,2,2)
heatmap(var_noise,N,log10(CRLB_a),'Colormap', parula(256),'ColorbarVisible','on','YLabel','N','XLabel','Noise Variance');
title('CRLB for a1 coefficient estimate')
savefig(figure(3),'figures/fig2_15.fig')
saveas(figure(3),'figures/forlatex/fig2_15','epsc')
%(??????????add colorbar label????????????)

%obtain a for CRLB estimates with an AR(1) model

%for standarised
[a_order1_std,e_order1_std] = aryule(stand_cp,1);

%for original data
[a_order1,e_order1] = aryule(closing_prices,1);

