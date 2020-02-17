%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2 Cross-correlation function%
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

N=1000;
x=randn(1,1000); %1000-sample WGN realisation
y=filter(ones(9,1),1,x);


%%  2.2.1 - unbiased estimate of CCF for x and y %% 
%------------------------------------------------------------------------------%
[ccf,timelag] = xcorr(y,x,'unbiased');

figure(1)
stem(timelag,ccf,'b','Linewidth',1.5);
title({'CCF for sequences x and y'},'Linewidth',1.5); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag (\tau)');    
xlim([-20 20]);
savefig(figure(1),'figures/fig2_6.fig')
saveas(figure(1),'figures/forlatex/fig2_6','epsc')
%%
%effect of filter order
clc
f_order=[5,15,30];
axis_lim=[20,20,40];
figure(2)
for subplot_index=1:length(f_order)
    subplot(1,length(f_order),subplot_index);
    
    y_2=filter(ones(f_order(subplot_index),1),1,x);
    
    ccf_y_2=xcorr(y_2, x,'unbiased');
    stem([-axis_lim(subplot_index):1:axis_lim(subplot_index)],ccf_y_2(-axis_lim(subplot_index)+1000:1:axis_lim(subplot_index)+1000),'r','.','Linewidth',1.1);
    str=sprintf('MA of order %d',f_order(subplot_index));title(str); 
    ylabel('Correlation, R(\tau)');
    xlabel('Time Lag (\tau)');
    xlim([-axis_lim(subplot_index) axis_lim(subplot_index)]);
end
savefig(figure(2),'figures/fig2_7.fig')
saveas(figure(2),'figures/forlatex/fig2_7','epsc')

%%

