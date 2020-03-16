%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 Cramer-Rao Lower Bound   %
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

% Testing the pgm function for WGN of different lengths

%For N values of [128,256,512]

%test wgn signals
% wgn_128 =randn(1,128);
% wgn_256 =randn(1,256);
% wgn_512 =randn(1,512);
load('wgn_randomsignals.mat') %for consistency

figure(1)
[pgm_128, freq_128]=pgm(wgn_128);
subplot(1,3,1)
plot(freq_128,pgm_128,'b','Linewidth',0.9);
xlabel('Normalised frequency');
ylabel('PSD estimate');
title('128-sample realisation of WGN');

[pgm_256, freq_256]=pgm(wgn_256);
subplot(1,3,2)
plot(freq_256,pgm_256,'b','Linewidth',0.9);
xlabel('Normalised frequency');
ylabel('PSD estimate');
title('256-sample realisation of WGN');

[pgm_512, freq_512]=pgm(wgn_512);
subplot(1,3,3)
plot(freq_512,pgm_512,'b','Linewidth',0.9);
xlabel('Normalised frequency');
ylabel('PSD estimate');
title('512-sample realisation of WGN');

% savefig(figure(1),'figures/fig3_1.fig')
% saveas(figure(1),'figures/forlatex/fig3_1','epsc')

%% --------------------------------------------------------------------------
% 3.1 - Averaged periodogram estimates

h=0.2*[1 1 1 1 1];% impulse response of the FIR filter
%filtering
smooth_pgm_128 = filter(h,1,pgm_128);
smooth_pgm_256 = filter(h,1,pgm_256);
smooth_pgm_512 = filter(h,1,pgm_512);

figure(2)
subplot(1,3,1)
plot(freq_128,smooth_pgm_128,'b');
xlabel('Normalised frequency');
ylabel('PSD estimate');
title({'Smoothened PSD of','128-sample realisation of WGN'});

subplot(1,3,2)
plot(freq_256,smooth_pgm_256,'b');
xlabel('Normalised frequency');
ylabel('PSD estimate');
title({'Smoothened PSD of','128-sample realisation of WGN'});

subplot(1,3,3)
plot(freq_512,smooth_pgm_512,'b');
xlabel('Normalised frequency');
ylabel('PSD estimate');
title({'Smoothened PSD of','128-sample realisation of WGN'});

% savefig(figure(2),'figures/fig3_2.fig')
% saveas(figure(2),'figures/forlatex/fig3_2','epsc')

%% -----------------------------------------------------------------------
%3.1.2 - PSD estimate for 8 non-overlapping 128-sample segments

%whole_wgn= randn(1,1024); %generate the WGN signal
load('wgn_1024.mat')
[whole_pgm, whole_freq]=pgm(whole_wgn); %periodogram of the whole
%variation of original series
mean_whole_pgm=mean(whole_pgm);
var_whole_pgm=var(whole_pgm);

figure(3)
for segment=1:8
    
    wgn(segment,:) = whole_wgn((((segment-1)*128)+1):(segment*128));
    
    %Computing the periodogram of the just created segment
    [pgm_wgn(segment,:), freq(segment,:)]=pgm(wgn(segment,:));
    
    subplot(2,4,segment)
    plot(freq(segment,:),pgm_wgn(segment,:),'b'); hold on;
    xlabel('Normalised frequency')
    ylabel('PSD estimate')
    str=sprintf('Segment: %d',segment);title(str);
    
    vars(segment)=var(pgm_wgn(segment,:));
    means(segment)=mean(pgm_wgn(segment,:));
    
end

% savefig(figure(3),'figures/fig3_3.fig')
% saveas(figure(3),'figures/forlatex/fig3_3','epsc')

%variation of psd estimate of series with segmentation
mean_segments=mean(means);
var_segments=mean(vars);

figure(5)
subplot(1,2,1)
stem(1:8,means,'Linewidth',2); hold on;
yline(mean_segments,'Linewidth',2);
yline(1,'r--','Linewidth',2);
title('Mean of each PSD estimate')
ylabel('PSD estimate');
xlabel('Segment index');xticks(1:8)
legend('\fontsize{11}Mean PSD of segment','\fontsize{11}Average of mean PSDs','\fontsize{11}Theoretical PSD')
legend('Location', 'southwest');

subplot(1,2,2)
stem(1:8,vars,'Linewidth',2); hold on
yline(var_segments,'Linewidth',2);
title('Variance of each PSD estimate')
ylabel('Variance');
xlabel('Segment index');xticks(1:8)
legend('\fontsize{11}Variance of PSD of segment','\fontsize{11}Average of PSD variances')
legend('Location', 'southwest');

% savefig(figure(5),'figures/fig3_5.fig')
% saveas(figure(5),'figures/forlatex/fig3_5','epsc')

%% -----------------------------------------------------------------------
%3.1.2 - PSD estimator: Averaged Periodogram (av_pmg)

av_pgm=mean(pgm_wgn,1); %average columns,i.e. accros segments

figure(4)
plot(whole_freq,whole_pgm,'y','Linewidth',1); hold on;
plot(freq_128,av_pgm,'b','Linewidth',2);
xlabel('Normalised frequency')
ylabel('PSD estimate')
title('Original and Averaged Periodogram')
set(gca,'Color','k')
legend('\fontsize{11}Original Periodogram','\fontsize{11}Averaged Periodogram')
legend('Color','w','Location','northeast')
grid on;


%variation of averaged periodogram
mean_av_pgm=mean(av_pgm);
var_av_pgm=var(av_pgm);














