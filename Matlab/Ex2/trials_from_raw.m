
%%  Obtain trials from ECG results  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

load 'ecg\RAWr'

trial_1=data(time>0.001563*(60*60*24) & time < 0.00281*(60*60*24));

trial_2=data(time>0.004303*(60*60*24) & time < 0.005575*(60*60*24));

trial_3=data(time>0.006927*(60*60*24) & time < 0.008316*(60*60*24));


%% convert to RRI
%[xRRI,fsRRI]=ECG_to_RRI(xECG,fsECG);

fsECG = 1000; %ecg sampling frequency


[xRRI_t1,fsRRI]=ECG_to_RRI(trial_1,fsECG);
% reject abnormalities after 1.231 and 10.115



[xRRI_t2,fsRRI]=ECG_to_RRI(trial_2,fsECG);
% reject abnormalities after 71.142 and 73.695



[xRRI_t3,fsRRI]=ECG_to_RRI(trial_3,fsECG);
% reject abnormalities after 1.122, 2.242, 3.745, 4.533, 5.338, 10.461,
% 26.957 and 43.123
