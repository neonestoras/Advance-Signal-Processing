
%%  1.1 Statistical estimation  %%
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

%%----------------------------------------------------
%create a signal x[n] of 1000 realisations
%that where each x[n] is a realisation of a 
%uniform random variable X ? U(0, 1) at time instant n

x=rand(1000,1);

figure(1)
plot(x);
title('Realisation of x[n] for 1000 samples');
xlabel('Sample number');ylabel('Amplitude');

%% Part1------------------------------------------------------
%compute the sample mean
sample_mean = mean(x);

%calculated mean of the uniform
%distribution from 0 to 1 is 0.5

%mean percentage error
mean_error_percentage=abs((sample_mean-0.5)/0.5)*100;

%% Part2----------------------------------------------------------
%Computing the sample standard deviation of x







