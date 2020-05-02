%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.2 Spectrum of autoregressive process   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Filter a 1064-sample WGN sequence x using the MATLAB command filter with
% the coef?cients b = 1 and a = [1, 0.9],

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
var(x)
mean(x)
var(y)
mean(y)

%% -----------------------------------------------------------------------%
%  3.2.1 - Cut off frequency of spectral estimate %% 

[h,w] = freqz(b, a, 512); % h:"EXACT psd"

% [H,W] = freqz(B,A,N) returns the N-point complex frequency response
%     vector H and the N-point frequency vector W in radians/sample of
%     the filter with numerator and denominator polynomial coefficients
%     in b and a, respectively.

%plotted in next part using plot(w/(2*pi),abs(h).^2)

% cu off frequency is where there is a drop of 3dB
exact_psd=abs(h).^2;
max_exact_psd=max(exact_psd);
Amplitude_at_3dBdrop=10^(-3/20)*max_exact_psd;

syms f
eqn= Amplitude_at_3dBdrop/max_exact_psd==1./((abs(1+0.9*exp(-2*pi*f))).^2);
f_coff=eval(abs(solve(eqn,f)));    

figure(5)
RGB2=[6 97 26];
f=0:0.01:1;
plot(f,20*log10(1./((abs(1+0.9*exp(-2*pi*f))).^2)),'Linewidth',2)
xline(f_coff,'--k','Linewidth',2)


%% -----------------------------------------------------------------------%
%  3.2.2 - Calculate the periodogram using pgm function and compare to the above%%
%  3.2.3 - 3.2.2 but zoomed in and ecplain the effect of rectangular window%%
RGB1=[0.0 0.4 0.0];
[pgm_y, freq_sun]=pgm(y);

%plot of Dirichlet function
figure (6)
samp_freq = 1000;
fr_ax = (-pi:2*pi/samp_freq:pi-2*pi/samp_freq);
for N=[5,10,20,50,100]
    clear Dirichlet
    Dirichlet=sinc(N.*fr_ax./2)./sinc(fr_ax./2);
    plot(fr_ax,(1/N).*Dirichlet.^2,'Linewidth',1.3);hold on;
end
title('$Nasinc_{N}^{2}(\omega)$','Fontsize',14,'Interpreter','latex','Linewidth',2)
legend('N=5','N=10','N=20','N=50','N=100','Fontsize',12);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (AU)')
xticks(-pi:pi/2:pi)
xticklabels({'- \pi',' - \pi / 2 ','0',' \pi / 2',' \pi ','Interpreter','latex'})
xlim([-pi,pi])
%  
% savefig(figure(6),'figures/fig3_dirichlet.fig')
% saveas(figure(6),'figures/forlatex/fig3_dirichlet','epsc')


%sinc function == DFT of rectangular window that affects the periodogram
figure(7)
fr_ax = (-pi:2*pi/samp_freq:pi-2*pi/samp_freq).*samp_freq./(2*pi);
for N=[5,10,20,50,100]
    plot(fr_ax,(1/N).*fftshift(abs(fft(rectwin(N),samp_freq))).^2,'Linewidth',1.5);
    hold on; 
end
title('$ \frac{1}{N} DFT \lbrace w_N[n] \rbrace^2 $','Fontsize',14,'Interpreter','latex','Linewidth',2)
xlabel('Frequency (rad)')
ylabel('Magnitude (AU)')
legend('N=5','N=10','N=20','N=50','N=100');
xlim([-200,200])

savefig(figure(7),'figures/fig3_sinc.fig')
saveas(figure(7),'figures/forlatex/fig3_sinc','epsc')


 %% -----------------------------------------------------------------------%
        % 3.2.4 - Calculate the model based PSD estimate

        [ACF_y,timelag] = xcorr(y,'unbiased');

        Ry0=ACF_y(timelag==0); %Ry(0)
        Ry1=ACF_y(timelag==1); %Ry(1)

        alpha1=-Ry1/Ry0;
        estimated_var=Ry0+alpha1*Ry1;

        % Calculate the model using the parameters found
        [h_model,w_model]=freqz(1,[1 alpha1],512);
        
        
        xlim([0 0.5]);

%% plot all 3.2.2, 3.2.3, 3.2.4

figure (2)
for subplot_index=1:2 
    subplot(1,2,subplot_index)
    plot(freq_sun,pgm_y,'r','Linewidth',1.2); hold on; %periodogram based psd
    xlabel('Normalised Frequency');
    ylabel('Magnitude');
    plot(w_model/(2*pi),estimated_var*(abs(h_model).^2),'b','Linewidth',2);%model based
    plot(w/(2*pi),abs(h).^2,'k','Linewidth',1.8);%exact
    
    if subplot_index==1
        title('PSD of filtered signal');
        xlim([0 0.5]);
    elseif subplot_index==2
        title('Zoomed PSD of filtered signal');
        xlim([0.4 0.5]);
    end
    xline(f_coff,'--','Linewidth',2)
    legend('Periodogram (using pgm)','Model based PSD','Exact PSD','Cut-off frequency','Location','northwest')
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


set(groot, 'defaultFigurePosition', [100, 100, 1200, 620]);
orders=[1 2 6 12 36]; %oders tested

for i=1:2 %repeat for standardised data
    %figure (2+i)
    figure(89)
    if i==1
        datasetlabel='original';
    elseif i==2
        datasetlabel='standardised';
    end
    for order_index=1:length(orders) 
        subplot(2,length(orders),order_index+length(orders)*(i-1))
        plot(freq_sun(i,:),pgm_sun(i,:));hold on;
        
         %to find the coefficient (a) and variance (e) required by the model
        [a,e] = aryule(suns(:,i), orders(order_index));
        [h_ar_model,w_ar_model]=freqz(1,a,512);
        plot(w_ar_model/(2*pi),e.*(abs(h_ar_model).^2),'Linewidth',1.5);
        xlabel('\fontsize{10}Normalised Frequency');
        ylabel('Magnitude');
        str=sprintf('AR(%d) of\n%s\nsunspot',orders(order_index),datasetlabel);
        title(str)
        
        xlim([0 0.5]); ylim([0 250000]);
        if order_index==1
            legend('\fontsize{9}Periodogram','\fontsize{9}Model-Based','Location','northeast');
        end
    end
end

% savefig(figure(89),'figures/fig3_89.fig')
% saveas(figure(89),'figures/forlatex/fig3_89','epsc')







