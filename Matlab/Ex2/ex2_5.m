%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.5 Realworldsignals: ECG from iAmp experiment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%  2.5.1 - Heart rate probability density estimate(PDE).  %% 
%-------------------------------------------------------------------------%

load 'xRRIr.mat'

%calculating the heart rate for each trial
hrate_1= 60./xRRI_t1;
hrate_2= 60./xRRI_t2;
hrate_3= 60./xRRI_t3;


%=======%
% Using trial 1 only, normal breathing
%average for every 10 samples (windows of 10) for smoother estimate
% using the selected coefficients a

a=[1 0.6];

avg_hr=[];

window=10; %window size
w=window/2; %for using in code

for i=1:((length(hrate_1)/w))
    avg_hr(i,1)=a(1)*mean(hrate_1(((i-1)*w+1):((i*w))));
    avg_hr(i,2)= a(2)*mean(hrate_1(((i-1)*w+1):((i*w))));
end

%Obtaining the PDE
% obtain probability density estimate F 
% using the normal kernel function

[F_1,XI_1]=ksdensity(avg_hr(:,1));
[F_2,XI_2]=ksdensity(avg_hr(:,2));
[F_3,XI_3]=ksdensity(hrate_1); 

figure (1)

area(XI_3,F_3,'LineWidth', 1.2,'FaceAlpha',0.7); hold on;
area(XI_1,F_1,'LineWidth', 1.2,'FaceAlpha',0.6);
area(XI_2,F_2,'LineWidth', 1.2,'FaceAlpha',0.75);
title('Probability Density Estimates of heart rates');
xlabel('Heart Rate (bpm)'); 
ylabel('Probability');
legend('Original HR','Average HR with a=1','Average HR with a=0.6')

% savefig(figure(1),'figures/fig2_16.fig')
% saveas(figure(1),'figures/forlatex/fig2_16','epsc')


%identify  peak values
peaks=zeros(2,3);
%1 for value , 2 for bpm
[peaks(1,1),I]=max(F_1);
peaks(2,1)=XI_1(I);
[peaks(1,2),I]=max(F_2);
peaks(2,2)=XI_2(I);
[peaks(1,3),I]=max(F_3);
peaks(2,3)=XI_3(I);

%and also std
hr_std=zeros(1,3);
hr_std(1)=std(avg_hr(:,1));
hr_std(2)=std(avg_hr(:,2));
hr_std(3)=std(hrate_1);

%%  2.5.2 - AR modelling of heart rate  %% 
%----------------------------------------%

%use detrend to ensure zero mean data

xRRI_t1_zm = detrend(xRRI_t1,'constant');
xRRI_t2_zm = detrend(xRRI_t2,'constant');
xRRI_t3_zm = detrend(xRRI_t3,'constant');

%compute the ACF of each
[ACF_t1,tlag_1] = xcorr(xRRI_t1_zm,'unbiased');
[ACF_t2,tlag_2] = xcorr(xRRI_t2_zm,'unbiased');
[ACF_t3,tlag_3] = xcorr(xRRI_t3_zm,'unbiased');

%colorsRGB
purplecol= [ 0.700 0.300 0.700];
orangecol=[ 0.9100 0.4100 0.1700];
%set(groot, 'defaultFigurePosition', [0, 100, 1600, 300]);

figure (2)
subplot(1,3,1)
plot(tlag_1,ACF_t1,'Linewidth',1.5);
title('ACF estimate of trial 1'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');
xlim([-200 200])

subplot(1,3,2)
plot(tlag_2,ACF_t2,'Color',orangecol,'Linewidth',1.5);
title('ACF estimate of trial 2'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');
xlim([-200 200])

subplot(1,3,3)
plot(tlag_3,ACF_t3,'Color',purplecol,'Linewidth',1.5);
title('ACF estimate of trial 3'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');
xlim([-200 200])

% savefig(figure(2),'figures/fig2_17.fig')
% saveas(figure(2),'figures/forlatex/fig2_17','epsc')

%compute the PACF of each
[PACF_t1,tlag_1] = parcorr(xRRI_t1_zm,'Method','Yule-Walker','NumLags',10);
[PACF_t2,tlag_2] = parcorr(xRRI_t2_zm,'Method','Yule-Walker','NumLags',10);
[PACF_t3,tlag_3] = parcorr(xRRI_t3_zm,'Method','Yule-Walker','NumLags',10);

figure (3)
subplot(1,3,1)
stem(tlag_1,PACF_t1,'Linewidth',1.9);hold on;
yline(1.96/sqrt(288), '--k','Linewidth',1.4); % 95% confidence interval
yline(-1.96/sqrt(288), '--k','Linewidth',1.4);
title('PACF estimate of trial 1'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');

subplot(1,3,2)
stem(tlag_2,PACF_t2,'Color',orangecol,'Linewidth',1.9);
yline(1.96/sqrt(288), '--k','Linewidth',1.4); % 95% confidence interval
yline(-1.96/sqrt(288), '--k','Linewidth',1.4);
title('PACF estimate of trial 2'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');

subplot(1,3,3)
stem(tlag_3,PACF_t3,'Color',purplecol,'Linewidth',1.8);
yline(1.96/sqrt(288), '--k','Linewidth',1.4); % 95% confidence interval
yline(-1.96/sqrt(288), '--k','Linewidth',1.4);
title('PACF estimate of trial 3'); 
ylabel('Correlation (R(\tau))');
xlabel('Time Lag (\tau)');

% savefig(figure(3),'figures/fig2_18.fig')
% saveas(figure(3),'figures/forlatex/fig2_18','epsc')


%identifying th model order using AIC MDL AICc criteria
mdl=[];
aic=[];
aicc=[];
N=zeros(1,3);
N(1)=length(xRRI_t1_zm);
N(2)=length(xRRI_t2_zm);
N(3)=length(xRRI_t3_zm);


for p= 1:10 %i=1 for a0, p: parameters iteration up to 10
    % trial 1
    [A,e(p,1)] = aryule(xRRI_t1_zm,p);
    mdl (:,p,1) = log10(e(p,1)) + ((p*log10(N(1)))/N(1)); % minimum description length criterion
    aic (:,p,1) = log10(e(p,1)) + (2*p)/N(1); %Akaike information criterion
    aic_c (:,p,1) =  aic(:,p,1) + (((2*p)*(p+1))/(N(1)-p-1)); %Corrected AIC 
    
     % trial 2
    [A,e(p,2)] = aryule(xRRI_t2_zm,p);
    mdl (:,p,2) = log10(e(p,2)) + ((p*log10(N(2)))/N(2)); % minimum description length criterion
    aic (:,p,2) = log10(e(p,2)) + (2*p)/N(2); %Akaike information criterion
    aic_c (:,p,2) =  aic(:,p,2) + (((2*p)*(p+1))/(N(2)-p-1)); %Corrected AIC 
    
     % trial 1
    [A,e(p,3)] = aryule(xRRI_t3_zm,p);
    mdl (:,p,3) = log10(e(p,3)) + ((p*log10(N(3)))/N(3)); % minimum description length criterion
    aic (:,p,3) = log10(e(p,3)) + (2*p)/N(3); %Akaike information criterion
    aic_c (:,p,3) =  aic(:,p,3) + (((2*p)*(p+1))/(N(3)-p-1)); %Corrected AIC 
end


figure (4)
for subplot_index=1:3
    subplot(1,3,subplot_index)
    plot(mdl(:,:,subplot_index),'b','LineWidth' ,2.2); hold on;
    plot(aic(:,:,subplot_index),'r','LineWidth' ,1.35);
    plot(aic_c(:,:,subplot_index),'LineWidth' ,1.5);
    plot(log10(e(:,subplot_index)),'k--','LineWidth' ,1.5);
    str=sprintf('MDL, AIC and AICc for Trial %d',subplot_index);
    title(str);
    str=sprintf('Model order\n');
    xlabel(str);xlim([1 10]);
    ylabel('Prediction Error');
end
legend('\fontsize{11}MDL', '\fontsize{11}AIC','\fontsize{11}AIC corrected','\fontsize{11}Loss function')
legend('Location','eastoutside','orientation','horizontal','Box','off');








