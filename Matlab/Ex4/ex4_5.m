%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.5 Speech Recognition     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% 4.5.1 - Voice signal predictor
[A,~] = audioread('Recorded_sounds\A.m4a'); 
[E,~] = audioread('Recorded_sounds\E.m4a');
[S,~] = audioread('Recorded_sounds\S.m4a');
[T,~] = audioread('Recorded_sounds\T.m4a');
[X,f_sampl] = audioread('Recorded_sounds\X.m4a');

N=1100; %length of saved sound
sounds=zeros(N,5);
sounds(:,1)=A(62500:62500+(N-1));
sounds(:,2)=E(53300:53300+(N-1));
sounds(:,3)=S(62000:62000+(N-1));
sounds(:,4)=T(75700:75700+(N-1));
sounds(:,5)=X(65000:65000+(N-1));
var_s=var(sounds,0,1);

letters=['A';'E';'S';'T';'X'];

figure(1)
mu_test=[0.01 0.02 0.04 0.06 0.12];
p_axis=1:40;
for l=1:length(letters)
    subplot(1,5,l)
    for m=1:length(mu_test)
        for p=1:40
            [~, e_48, w_48] = adaptive_lms(sounds(:,l),mu_test(m),p);
            
            cse(p,m)=sum(e_48.^2); %cummulative squared error
            
            R(p,m)=10*log10(var_s(l)/var(e_48)); %predictor gain
            
            
            
        end
    end
    
    for p=1:40
        [~, e_gs, w_gs, ~] = adaptive_gs_lms(sounds(:,l), 0.1, 0.125, p);
        cse_gs(p)=sum(e_gs.^2);
        R_gs(p)=10*log10(var_s(l)/var(e_gs)); 
    end
    
    plot(p_axis,20*log10(cse),'Linewidth',1.35);hold on;
    plot(p_axis,20*log10(cse_gs),'k','Linewidth',1.4);
    title(['Sound:  ' letters(l)]);
    ylabel('CSE (dB)');
    xlabel('Model Order');
    xlim([1 40])
end
legend('\mu=0.01', '\mu=0.02','\mu=0.04','\mu=0.06','\mu=0.12','GS(\mu=0.1)');
legend('Fontsize',11)

% savefig(figure(1),'figures/fig4_7.fig')
% saveas(figure(1),'figures/forlatex/fig4_5','epsc')

Prediction_gain_48kHz=table({'order','\mu=0.01', '\mu=0.02','\mu=0.04','\mu=0.06','\mu=0.12','GS(\mu=0.1)'}',[1:40; R';R_gs]);






%% -----------------------------------------------------------------------%
% 4.5.3 - change sampling freq
downsampling_factor=3;
f_sampl=f_sampl/downsampling_factor;
A_d=downsample(A,3);
E_d=downsample(E,3);
S_d=downsample(S,3);
T_d=downsample(T,3);
X_d=downsample(X,3);


N=1000;
d_sounds=zeros(N,5);
d_sounds(:,1)=A_d(floor(62500/3):floor(62500/3)+(N-1));
d_sounds(:,2)=E_d(floor(53300/3):floor(53300/3)+(N-1));
d_sounds(:,3)=S_d(floor(62000/3):floor(62000/3)+(N-1));
d_sounds(:,4)=T_d(floor(75700/3):floor(75700/3)+(N-1));
d_sounds(:,5)=X_d(floor(65000/3):floor(65000/3)+(N-1));
var_d_s=var(d_sounds,0,1);


figure(2)
mu=0.1;
for l=1:length(letters)
    subplot(1,5,l)
    
    for p=1:40
            [~, e_48, ~] = adaptive_lms(sounds(:,l),mu,p);
            [~, e_16, ~] = adaptive_lms(d_sounds(:,l),mu,p);
            
            cse(p,1)=sum(e_48.^2);
            cse(p,2)=sum(e_16.^2);%cummulative squared error
            R_48(p,l)=10*log10(var_s(l)/var(e_48));
            R_16(p,l)=10*log10(var_d_s(l)/var(e_16));%predictor gain
    end
    plot(p_axis,R_48(:,l),'Linewidth',1.35);hold on;
    plot(p_axis,R_16(:,l),'Linewidth',1.35);
    title(['Sound:  ' letters(l)]);
    ylabel('Prediction Gain');
    xlabel('Model Order');
    xlim([1 40])
end
legend('f_s=48kHz','f_s=16kHz');
legend('Fontsize',12)

% savefig(figure(2),'figures/fig4_8.fig')
% Prediction_gain_sampling_freq=table({'sound','f_s=48kHz','f_s=16kHz'},[letters'; num2str(max(R_48)'); num2str(max(R_48)')]);













