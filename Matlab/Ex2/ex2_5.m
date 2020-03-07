%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.5 Realworldsignals: ECG from iAmp experiment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%  2.5.1 - Heart rate probability density estimate(PDE).  %% 
%-------------------------------------------------------------------------%

load 'xRRIr.mat'

%calculating the heart rate for each trial
hrate_1= 60./xRRI_t1;
hrate_2= 60./xRRI_t2;
hrate_3= 60./xRRI_t3;


%=======%
% Using trial 1 only, normal breathing
%average for every 10 samples (windows of 10) for smoother estimate
% using the selected coefficients a

a=[1 0.6];

avg_hr=[];

window=10; %window size
w=window/2; %for using in code

for i=1:((length(hrate_1)/window))
    avg_hr(i,1)=a(1)*mean(hrate_1(((i-1)*w+1):((i*w))));
    avg_hr(i,2)= a(2)*mean(hrate_1(((i-1)*w+1):((i*w))));
end

%Obtaining the PDE
[f_a1,xi_a1]=ksdensity(avg_hr(1));
[f_a2,xi_a2]=ksdensity(avg_hr(2));



