%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.5 Realworldsignals: Respiratory sinus arrhythmia from RR-Intervals   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 400]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);
RGB1=[1, 70, 1]/100;

%Load RRI data used in ex2.5
load 'xRRIr.mat'
%remove their offset
dRRI_1=detrend(xRRI_t1);
dRRI_2=detrend(xRRI_t2);
dRRI_3=detrend(xRRI_t3);


%% -----------------------------------------------------------------------%
% 3.5.a - PSD estimaation of RRI

%cargulate periodograms
[pgm_t1, freq_t1]=pgm(dRRI_1);
[pgm_t2, freq_t2]=pgm(dRRI_2);
[pgm_t3, freq_t3]=pgm(dRRI_3);

figure(1)
for subplot_index=1:3
    subplot(3,3,subplot_index)
    if subplot_index==1
        plot(freq_t1,pgm_t1,'g')
    elseif subplot_index==2
        plot(freq_t2,pgm_t2,'g')
    elseif subplot_index==3
        plot(freq_t3,pgm_t3,'g')
    end
    xlabel('Normalised Frequency')
    ylabel('PSD')
    xlim([0 0.15])
    title(['Periodogram of Trial ',num2str(subplot_index)])
    set(gca,'Color','k')
end

%% Averaged periodogram
segments=[2 6];
for subplot_index = 1:length(segments)
    
    num_segments=segments(subplot_index);
    L_segment=floor(length(dRRI_1)/num_segments);
    segm_1=[];
    segm_2=[];
    segm_3=[];
    for s=1:num_segments
        segm_1(s,:) = dRRI_1((((s-1)*L_segment)+1):(s*L_segment)); %devide to segments
        segm_2(s,:) = dRRI_2((((s-1)*L_segment)+1):(s*L_segment));
        segm_3(s,:) = dRRI_3((((s-1)*L_segment)+1):(s*L_segment));
        [psd_segm_1(s,:),f_segm_1(s,:)]=pgm(segm_1(s,:)); %periodogram of segment
        [psd_segm_2(s,:),f_segm_2(s,:)]=pgm(segm_2(s,:));
        [psd_segm_3(s,:),f_segm_3(s,:)]=pgm(segm_3(s,:));
    end
%     Averaged periodogram
%     This will take the average value of each column of our matrix
    avg_pgm(1,:) = mean(psd_segm_1,1);
    avg_pgm(2,:) = mean(psd_segm_2,1);
    avg_pgm(3,:) = mean(psd_segm_3,1);
    
    for t=1:3
        subplot(3,3,(subplot_index)*3+t)
        plot(f_segm_1(1,:),avg_pgm(t,:),'y')
        xlabel('Normalised Frequency')
        ylabel('PSD')
        str=sprintf('Window Length %d\n Trial %d  ',L_segment,t);
        title(str)
        xlim([0 0.15])
        set(gca,'Color','k')
    end

    clear segm_1
    clear segm_2
    clear segm_3
    clear psd_segm_1
    clear psd_segm_2
    clear psd_segm_3
    clear f_segm_1
    clear f_segm_2
    clear f_segm_3
    clear avg_pgm
end


