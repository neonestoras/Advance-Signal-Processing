%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3 Autoregressive modelling %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%  2.3.1 - AR(2) coefficients %% 
%-------------------------------------------------------------------------%
% Generate 100 samples of the uniformly distributed random variables
% a1 in [-2:5; 2:5]
% and a2 in [-1:5; 1:5].
% generate N random numbers in the interval (a,b) with the formula 
%r = a + (b-a).*rand(N,1)

L=1000;%process length
N_tested=[100 5000];

for N_index=1:length(N_tested)

N=N_tested(N_index);%number of samples
a=zeros(2,N);

figure (1)
subplot(1,length(N_tested),N_index)  

a(1,:) = -2.5 + (2.5+2.5)*rand(N,1);
a(2,:) = -1.5 + (1.5+1.5)*rand(N,1);
wgn = randn(L,1); %white gaussian noise


%generate the signal x[n] with the give equation
x=zeros(N,L);
for i=1:N
   x(i,:) = filter(1,[1,-a(1,i),-a(2,i)],wgn);%minus sign for AR part coeff
end
%Obtaining the maximum values of each realisation of signal.
max_values = max(x,[],2);
%after observation of max_values the decided empirical bound is 1000
bound=1000;%anything above 3 figures

un_index=1;
stable_index=1;

for i=1:N
    if abs(x(i,:))<=bound
       stable_pairs(stable_index,1)= a(1,i);
       stable_pairs(stable_index,2)= a(2,i);
       stable_index=stable_index+1;
    else
        un_pairs(un_index,1)= a(1,i);
        un_pairs(un_index,2)= a(2,i);
        un_index=un_index+1;
    end
end



%stable pairs
plot(stable_pairs(:,1),stable_pairs(:,2),'*'); hold on;

%unstable pairs
plot(un_pairs(:,1),un_pairs(:,2),'x r')
str=sprintf('%d Pairs',N);title(str,'Linewidth',1.5); 
ylabel('a2');
xlabel('a1');


axis([-2.5 2.5 -1.5 1.5]);


%Plotting the Stability Triangle
%and the line separatig real and complex roots
lin=[-2:0.01:2];
plot(lin,-((lin).^2)/4,'--k','Linewidth',1.3);
plot([-2 2],[-1 -1],'k','Linewidth',1.3);
plot([-2 0],[-1 1],'k','Linewidth',1.3);
plot([0 2],[1 -1],'k','Linewidth',1.3);

if N==N_tested(length(N_tested))
    legend('\fontsize{10}Stable pairs','\fontsize{10}Unstable Pairs','\fontsize{10}Real-Imag Roots line','\fontsize{10}Stability Triangle')
    legend('Location','eastoutside');
else
    clear stable_pairs;
    clear un_pairs;
end

end
% savefig(figure(1),'figures/fig2_8.fig')
% saveas(figure(1),'figures/forlatex/fig2_8','epsc')


%%  2.3.2 - AR of sunspot, real world data %% 
%-------------------------------------------------------------------------%

load sunspot.dat;
%This is a 288x2 matrix (sunspot) from matlab
%column 1 is the year index and the second is the data

datalengths=[5 20 250];

sunspot_5= sunspot (1:datalengths(1),2);
[acf_sun_5,timelag_5]= xcorr(sunspot_5,'unbiased');
mean_sunspot_5=mean(sunspot_5);
sunspot_zm_5= sunspot_5 - mean_sunspot_5;
[acf_sun_zm_5,timelag_zm_5]= xcorr(sunspot_zm_5,'unbiased');


sunspot_20 = sunspot (1:datalengths(2),2);
[acf_sun_20,timelag_20] = xcorr(sunspot_20,'unbiased');
mean_sunspot_20 = mean(sunspot_20);
sunspot_zm_20= sunspot_20 - mean_sunspot_20;
[acf_sun_zm_20,timelag_zm_20]= xcorr(sunspot_zm_20,'unbiased');

sunspot_250 = sunspot (1:datalengths(3),2);
[acf_sun_250,timelag_250] = xcorr(sunspot_250,'unbiased');
mean_sunspot_250 = mean(sunspot_250);
sunspot_zm_250= sunspot_250 - mean_sunspot_250;
[acf_sun_zm_250,timelag_zm_250]= xcorr(sunspot_zm_250,'unbiased');

figure(2)
subplot(2,2,1)
stem(timelag_5,acf_sun_5,'.','b','Linewidth',1.2); hold on;
stem(timelag_zm_5,acf_sun_zm_5,'.','r','Linewidth',1.2)
title('ACF using 5 samples'); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag, \tau');
xlim([-4 4]);
hold off

subplot(2,2,2)
stem(timelag_20,acf_sun_20,'.','b','Linewidth',1.2); hold on;
stem(timelag_zm_20,acf_sun_zm_20,'.','r','Linewidth',1.2)
title('ACF using 20 samples sunspot'); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag, \tau'); 
xlim([-19 19])

subplot(2,2,3)
stem(timelag_250,acf_sun_250,'.','b','Linewidth',1.2); hold on;
stem(timelag_zm_250,acf_sun_zm_250,'.','r','Linewidth',1.2)
title('ACF using 250 samples'); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag, \tau');
xlim([-249 249])

subplot(2,2,4)
stem(timelag_250,acf_sun_250,'.','b','Linewidth',1.2); hold on;
stem(timelag_zm_250,acf_sun_zm_250,'.','r','Linewidth',1.2)
title('ACF using 250 samples, zoomed in'); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag, \tau');
xlim([-50 50])
legend('\fontsize{10}Original','\fontsize{10}Zero-mean','Location','South')
legend('Orientation','horizontal')




%%  2.3.3 - Yule Walker equations %% 
%-------------------------------------------------------------------------%

%First for the original sunspot series

%Calculate all reflection coefficients up to p=10


for i =1:10
    [a_or,e_or(i)] = aryule(sunspot(:,2),i); 
    display(a_or)
end
%also calculating the reflection coeff (rc)
[a_or,e_or,rc_or] = aryule(sunspot(:,2),10);
pacf_sun_org=-rc_or; %partiaal acf equal to rc*-1

%Secondly the standardised sunspot series (zero mean and unit variance)
stand_sunspot = (sunspot(:,2) - mean(sunspot(:,2)))./std(sunspot(:,2));
[a_stand,e_stand,rc_stand] = aryule(stand_sunspot,10);
pacf_sun_stand=-rc_stand;

figure(3)

stem(pacf_sun_org,'b','Linewidth',1.4); hold on;
stem(pacf_sun_stand,'r', 'MarkerFaceColor','red','MarkerEdgeColor','none','Linewidth',1.3);
yline(1.96/sqrt(288), '--g','Linewidth',1.1);  % 95% confidence interval 
yline(2.575/sqrt(288), '--m','Linewidth',1.1); % 99% confidence interval
yline(3/sqrt(288),'--k','Linewidth',1.1); % 3 std errors


yline(-1.96/sqrt(288), '--g','Linewidth',1.1);
yline(-2.575/sqrt(288), '--m','Linewidth',1.1);
yline(-3/sqrt(288),'--k','Linewidth',1.1);

title('PACF of the original and standarised sunspot series model'); 
ylabel('Correlation');
xlabel('Time-Lag (k)');
legend('\fontsize{10}PACF original series','\fontsize{10}PACF standarised series','\fontsize{10}95% confidence interval','\fontsize{10}99% confidence interval','\fontsize{10}3 standard errors');
legend('Location','eastoutside');

%desplay the order 2 coefficients of the standarised series
[a_stand_2,e_stand]=aryule(stand_sunspot,2);
display(a_stand_2)
%%  2.3.4 - Find the correct order of model for the standardised data %% 
%-------------------------------------------------------------------------%

N=length(stand_sunspot);

mdl=[];
aic=[];
aic_c=[];

for p= 1:10 %i=1 for a0, p: parameters iteration
    [A,e(p)] = aryule(stand_sunspot,p);
    mdl (:,p) = log10(e(p)) + ((p*log10(N))/N); % minimum description length criterion
    aic (:,p) = log10(e(p)) + (2*p)/N; %Akaike information criterion
    aic_c (:,p) =  aic(:,p) + (((2*p)*(p+1))/(N-p-1)); %Corrected AIC 
end

figure(4)
plot(mdl,'b','LineWidth' ,2.2); hold on;
plot(aic,'r','LineWidth' ,1.35);
plot(aic_c,'LineWidth' ,1.5);
plot(log10(e),'k--','LineWidth' ,1.5)
legend('MDL', 'AIC','AIC corrected','Loss function','Location','eastoutside')
title('MDL, AIC and AICc for standarised sunspot series')
xlabel('Model order');xlim([1 10])
ylabel('Prediction Error');

% savefig(figure(4),'figures/fig2_11.fig')
% saveas(figure(4),'figures/forlatex/fig2_11','epsc')

% The correct oder of the model is at the lowest points on the graph
%the !index! of the minimum point is the right order

[optimal_mdl,min_mdl_index]=min(mdl);
[optimal_aic,min_aic_index]=min(aic);
[optimal_aicc,min_aicc_index]=min(aic_c);
min_mdl_index
min_aic_index
min_aicc_index

%%  2.3.5 - Predicting using an AR model %% 
%-------------------------------------------------------------------------%
orders=[1 2 10];
predict_hor= [1 2 5 10];
predictions=[];
for subplot_index= 1:length(orders)
    model= ar(stand_sunspot,orders(subplot_index),'yw');
    %approach -> 'yw' : The Yule-Walker method
    for i= 1:length(predict_hor)
        predictions(i+4*(subplot_index-1),:)=predict(model,stand_sunspot,predict_hor(i));
    end
    
end

figure(5)

for subplot_index= 1:length(predict_hor)
    
    subplot(length(predict_hor),1,subplot_index)
    plot(stand_sunspot,'Linewidth',1.8); hold on;
    plot(predictions(1+(subplot_index-1),:),'Linewidth',1.2);
    plot(predictions(5+(subplot_index-1),:),'Linewidth',1.2);
    plot(predictions(9+(subplot_index-1),:),'Linewidth',1.2);
    str=sprintf('\n\n Prediction Horison = %d',predict_hor(subplot_index));
    title(str)
    xlabel('Sample'); ylabel('Magnitude')
    
end
legend('Original','AR(1)','AR(2)', 'AR(10)');
legend('Location','northoutside','Orientation','horizontal')

% set(groot, 'defaultFigurePosition', [100, 100, 1500, 460]);

figure(6)
i=1;
for order_index=1:length(orders)
    for hor_index= 1:length(predict_hor)
        subplot(3,4,4*(order_index-1)+hor_index)
    plot(stand_sunspot,'k','Linewidth',1.3); hold on;
    plot(predictions(i,:),'r','Linewidth',1.1);
    i=i+1;
    str=sprintf('Moder AR(%d)\nPrediction Horison = %d',orders(order_index),predict_hor(hor_index));
    title(str,'FontSize', 11)
    xlabel('Sample'); ylabel('Magnitude')
    end
end
legend('\fontsize{10}Standarised Original Model','\fontsize{10}Standarised Predicted Model')
legend('Location','southoutside','orientation','horizontal')
legend('Box','off')

%%prediction errors
mse_1=zeros(1,4);
mse_2=zeros(1,4);
mse_10=zeros(1,4);
for i=1:length(predict_hor)*length(orders)
    
    mse= immse(stand_sunspot,predictions(i,:)');
    
        if i<5
            mse_1(i) = mse;
    
        elseif (i>4) && (i<9)
            mse_2(i-4) = mse;
        
        elseif (i>8) && (i<13)
            mse_10(i-8) = mse;
        
        end
    
end

figure(7)
bar([mse_1; mse_2; mse_10]');hold on;
title('Mean Square Error for models of different orders');
xlabel('Prediction Horizon (years)')
xticklabels({'1','2','5','10'})
ylabel('MSE')
legend('order 1','order 2','order 10','Location','eastoutside');

% savefig(figure(7),'figures/fig2_13.fig')
% saveas(figure(7),'figures/forlatex/fig2_13','epsc')













