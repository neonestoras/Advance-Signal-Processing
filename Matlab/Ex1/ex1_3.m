%%  1.3 Estimation of probability distributions  %%
%-------------------------------------------------%
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


%% Part1------------------------------------------------------

% Test the pdf estimator function for a stationary process 
% using Gaussian pdf with different data length(N)

%create a histogram for comparison
N=1000;
v=randn(1,N);
figure (1);
hist(v,100);


figure (2);
for subplot_index=1:4
    subplot(2,2,subplot_index)
    N=10^(subplot_index+1);
    v=randn(1,N);
    
    pdf(v);hold on;
    xlabel('v values');xlim([-4,4]); ylabel('Estimate of pdf of v');
    str = sprintf('pdf estimation using\ndata length of %d and 100 Bins',N);
    title(str);
    
    axis=[-4:0.001:4];
    plot(axis,normpdf(axis,0,1),'r','Linewidth',1.5);
    legend('\fontsize{10}Estimated pdf','\fontsize{10}Actual pdf','Location','northwest');
    hold off;
end

%% Part2------------------------------------------------------

%Approximatethe pdf for ergodic and stationary processes studied in Part 1.2 
%with data lengths of 100, 1000 and 10000]
%Compare estimated densities

% !Only preocess rp3 is both ergotic and stationary!

M=1;%single realisation used
figure (3);
subplot_index=1;
for N=[100,1000,10000,100000]
    
    subplot(1,4,subplot_index)
    
    v2=rp3(M,N);
    pdf(v2);hold on;
    xlabel('sample value'); ylabel('pdf estimate');
    str=sprintf('Data length %d',N);title(str);
    %Plotting the theoretical pdf
    line([-1,2],[(1/3),(1/3)],'Color','red','Linewidth',1.2);
    legend('\fontsize{10}Estimated pdf','\fontsize{10}Actual pdf','Location','southwest');
    
    subplot_index=subplot_index+1;
end

%% Part3------------------------------------------------------









