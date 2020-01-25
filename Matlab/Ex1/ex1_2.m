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

M1=100; %members of enseble
N1=99; %length of each
stds=zeros(3,M1); 
means=zeros(3,M1);
dtime_ax = [0:0.001:N1]; %discrete time axis
outputs=zeros(M1,N1,3);
outputs(:,:,1)= rp1(M1,N1);
outputs(:,:,2)= rp2(M1,N1);
outputs(:,:,3)= rp3(M1,N1);

figure(1)
for subplot_index=1:3
    subplot(1,3,subplot_index)
    means(subplot_index,:)=mean(outputs(:,:,subplot_index),2);
    plot(means(1,subplot_index),'Linewidth',1);hold on;
    %Plotting the theoretical mean
    if subplot_index==1
        plot(dtime_ax,0.02*dtime_ax,'r','Linewidth',1);
    elseif subplot_index==2
        line([0,100],[0.5,0.5],'Color','red','Linewidth',1);
    elseif subplot_index==3
        line([0,100],[0.5,0.5],'Color','red','Linewidth',1);
    end
    xlabel('Discrete time'),  ylabel('\fontsize{14}Magnitude')
    str = sprintf('Ensemble mean vs Time\nusing Stochastic Process %d', subplot_index);
    title(str);hold off;
end

























