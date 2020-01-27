%%  1.2 Stochastic processes  %%
%--------------------------------%
clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 300]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%-------------------------------------------------------

% explain the differences between the time averages and 
% ensemble averages, together with
% the stationarity and ergodicity of the process generated




%% Part1------------------------------------------------------
% Compute and plot the ensemble mean and standard deviation for each process
% Using M = 100 members of the ensemble,each of length N = 100

M1=100; %members of enseble
N1=100; %length of each
stds=zeros(3,N1); 
means=zeros(3,N1);
dtime_ax = [1:1:N1]; %discrete time axis
outputs=zeros(M1,N1,3);
outputs(:,:,1)= rp1(M1,N1);
outputs(:,:,2)= rp2(M1,N1);
outputs(:,:,3)= rp3(M1,N1);

figure(1)
for subplot_index=1:3
    means(subplot_index,:)=mean(outputs(:,:,subplot_index),1);
    subplot(2,3,subplot_index)
    plot(dtime_ax,means(subplot_index,:),'b','Linewidth',1);hold on;
    %Plotting the theoretical mean
    if subplot_index==1
        plot(dtime_ax,0.02*dtime_ax,'r','Linewidth',0.8);
    elseif subplot_index==2
        line([0,100],[0.5,0.5],'Color','red','Linewidth',0.8);
    elseif subplot_index==3
        line([0,100],[0.5,0.5],'Color','red','Linewidth',0.8);
    end
    xlabel('Discrete time'),  ylabel('\fontsize{14}Magnitude')
    str = sprintf('Ensemble mean vs Time\nusing Stochastic Process %d', subplot_index);
    title(str);hold off;
    
    stds(subplot_index,:)=std(outputs(:,:,subplot_index),1);
    subplot(2,3,subplot_index+3)
    plot(dtime_ax,stds(subplot_index,:),'b','Linewidth',1);hold on;
    %Plotting the theoretical std
    if subplot_index==1
        plot(dtime_ax,(5/sqrt(12))*sin((dtime_ax*pi)/N1),'r','Linewidth',0.8); 
    elseif subplot_index==2
        line([0,100],[(1/3),(1/3)],'Color','red','Linewidth',0.8);
    elseif subplot_index==3
        line([0,100],[(sqrt(3)/2),(sqrt(3)/2)],'Color','red','Linewidth',0.8);
    end
    xlabel('Discrete time'), ylabel('Magnitude')
    str = sprintf('Ensemble St.D. vs Time\nusing Stochastic Process %d', subplot_index);
    title(str);hold off;
    
end
saveas(figure(1),'figures/fig1_2_1','epsc');
%% Part2 ---------------------------------------------------------------------
% Generate M = 4 realisations of length N = 1000, 
% calculate the mean and standard deviation for each realisation.
M2=4;
N2=1000;
stds2=zeros(3,M2); 
means2=zeros(3,M2);

outputs2=zeros(M2,N2,3);
outputs2(:,:,1)= rp1(M2,N2);
outputs2(:,:,2)= rp2(M2,N2);
outputs2(:,:,3)= rp3(M2,N2);

for i=1:1:3
means2(i,:)=mean(outputs2(:,:,i)');
stds2(i,:)=std(outputs2(:,:,i)');
end

P1=table([1:1:M2]',means2(1,:)',stds2(1,:)');
P2=table([1:1:M2]',means2(2,:)',stds2(2,:)');
P3=table([1:1:M2]',means2(3,:)',stds2(3,:)');


T=table([1:1:M2]',means2(1,:)',stds2(1,:)',means2(2,:)',stds2(2,:)',means2(3,:)',stds2(3,:)')









