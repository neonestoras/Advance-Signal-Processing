%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.2 Spectrum of autoregressive process   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 350]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);



%% -----------------------------------------------------------------------%

%x=randn(1,1064); %1064-sampe WGN 
load('wgn_1064.mat')
b=1;
a=[1 0.9];
y=filter(b,a,x); 
y=y(41:length(y)); % due to transient filter effects first 40 samples are removed

figure (1)
subplot(2,1,1)
plot(x,'b')
xlabel('Sample');xlim([0 length(x)]);
ylabel('Magnitude');
title('Input Signal')

subplot(2,1,2)
plot(y,'r')
xlabel('Sample');xlim([0 length(y)]);
ylabel('Magnitude');
title('Filtered, output Signal')

% savefig(figure(1),'figures/fig3_6.fig')
% saveas(figure(1),'figures/forlatex/fig3_6','epsc')

%% -----------------------------------------------------------------------%
%  3.2.1 - Cut off frequency of spectral estimate %% 

[h,w] = freqz(b, a, 512);

% [H,W] = freqz(B,A,N) returns the N-point complex frequency response
%     vector H and the N-point frequency vector W in radians/sample of
%     the filter with numerator and denominator polynomial coefficients
%     in b and a, respectively.

%plotted tin next part using plot(w/(2*pi),abs(h).^2)

%% -----------------------------------------------------------------------%
%  3.2.2 - Calculate the periodogram using pgm function and compare to the above%%
%  3.2.3 - 3.2.2 but zoomed in%%
RGB1=[0.0 0.4 0.0];
[pgm_y, freq_sun]=pgm(y);

figure (2)
for subplot_index=1:2 
    subplot(1,2,subplot_index)
    plot(freq_sun,pgm_y,'r','Linewidth',1.2); hold on;
    xlabel('Normalised Frequency');
    ylabel('Magnitude');
    if subplot_index==1
        title('PSD of filtered signal');
        


        %% -----------------------------------------------------------------------%
        % 3.2.4 - Calculate the model based PSD estimate

        [ACF_y,timelag] = xcorr(y,'unbiased');

        Ry0=ACF_y(timelag==0); %Ry(0)
        Ry1=ACF_y(timelag==1); %Ry(1)

        alpha1=-Ry1/Ry0;
        estimated_var=Ry0+alpha1*Ry1;

        % Calculate the model using the parameters found
        [h_model,w_model]=freqz(1,[1 alpha1],512);
        
        plot(w_model/(2*pi),estimated_var*(abs(h_model).^2),'b','Linewidth',2);
        xlim([0 0.5]);
    elseif subplot_index==2
        title('Zoomed PSD of filtered signal');
        plot(w_model/(2*pi),estimated_var*(abs(h_model).^2),'b','Linewidth',2);
        xlim([0.4 0.5]);
    end
    plot(w/(2*pi),abs(h).^2,'k','Linewidth',1.8);
    legend('Periodogram (using pgm)','Model based PSD','PSD of AR(1) process','Location','northwest')
end


% savefig(figure(2),'figures/fig3_7.fig')
% saveas(figure(2),'figures/forlatex/fig3_7','epsc')

%% -----------------------------------------------------------------------%
% 3.2.4 - Repeat for Sunspot time Series

load sunspot.dat; %Data obtained from MATLAB itself
                    %(used in previous parts)

suns=zeros(288,2);
suns(:,1)=sunspot(:,2); %original data
suns(:,2)=detrend(sunspot(:,2)); %standardised data

% Estimate Periodogram using pmg function
pgm_sun=[]; freq_sun=[];
[pgm_sun(1,:), freq_sun(1,:)]=pgm(suns(:,1)');
[pgm_sun(2,:), freq_sun(2,:)]=pgm(suns(:,2)');


set(groot, 'defaultFigurePosition', [100, 100, 1200, 300]);
orders=[1 2 6 12 24];

for i=1:2 %repeat for standardised data
    figure (2+i)
    if i==1
        datasetlabel='original';
    elseif i==2
        datasetlabel='standardised';
    end
    for order_index=1:length(orders) 
        subplot(1,length(orders),order_index)
        plot(freq_sun(i,:),pgm_sun(i,:));hold on;
        
         %to find the coefficient (a) and variance (e) required by the model
        [a,e] = aryule(suns(:,1), orders(order_index));
        [h_ar_model,w_ar_model]=freqz(1,a,512);
        plot(w_ar_model/(2*pi),e.*(abs(h_ar_model).^2),'Linewidth',1.5);
        xlabel('\fontsize{10}Normalised Frequency');
        ylabel('Magnitude');
        str=sprintf('AR(%d) of\n%s\nsunspot',orders(order_index),datasetlabel);
        title(str)
        
        xlim([0 0.5]); ylim([0 250000]);
        legend('\fontsize{9}Periodogram','\fontsize{9}Model-Based','Location','northeast')
    end
end

% savefig(figure(3),'figures/fig3_8.fig')
% saveas(figure(3),'figures/forlatex/fig3_8','epsc')
% 
% savefig(figure(4),'figures/fig3_9.fig')
% saveas(figure(4),'figures/forlatex/fig3_9','epsc')








