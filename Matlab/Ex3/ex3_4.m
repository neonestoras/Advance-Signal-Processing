%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.4 Spectrogramfortime-frequencyanalysis: dialtonepa   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 310]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);
RGB1=[1, 70, 1]/100;


load('dialpad_frequencies') %first roow is the ASCII table number of the dialed character
                           % 2nd and 3rd are the frequencies


%% -----------------------------------------------------------------------%
% 3.4.1 - random London landline number

% line_num=[0 2 0 randi([0,9],1,8)]; %generate random number
load('random_landline_number')
sampl_f=32768;
segment_duration=0.25;
Lseg=sampl_f*segment_duration;
y= zeros(1,21*Lseg); % discrete time sequence



n=1:Lseg;
d_sequence=zeros(length(line_num)+1,Lseg);
idle_index=length(line_num)+1;
d_sequence(idle_index,:)=zeros(1,Lseg);  % idle sequence
for d=1:length(line_num)
    digit=char(line_num(1,d));
    f1=dialpad(2,(dialpad(1,:)==double(digit)));
    f2=dialpad(3,(dialpad(1,:)==double(digit)));
   
    d_sequence(d,:)=sin(2*pi*f1*(n-1)*(1/sampl_f)) + sin(2*pi*f2*(n-1)*(1/sampl_f)); 
end

y=[d_sequence(1,:) d_sequence(idle_index,:) d_sequence(2,:) d_sequence(idle_index,:) d_sequence(3,:) d_sequence(idle_index,:) d_sequence(4,:) d_sequence(idle_index,:) d_sequence(5,:) d_sequence(idle_index,:) d_sequence(6,:) d_sequence(idle_index,:) d_sequence(7,:) d_sequence(idle_index,:) d_sequence(8,:) d_sequence(idle_index,:) d_sequence(9,:) d_sequence(idle_index,:) d_sequence(10,:) d_sequence(idle_index,:) d_sequence(11,:)];
t_axis=0:1/sampl_f:(21*0.25)-(1/sampl_f);
set(groot, 'defaultFigurePosition', [100, 100, 850, 300]);
figure(1)
plot(t_axis,y,'Color',RGB1,'Linewidth',0.2);
title({'Discrete Time Sequence of  Landline number'; num2str(line_num)});
xlabel('Time (s)'); xlim([0,21*0.25])
ylabel('Amplitude')

% sound(y,sampl_f)
% savefig(figure(1),'figures/fig3_14.fig')
% saveas(figure(1),'figures/forlatex/fig3_14','epsc')





set(groot, 'defaultFigurePosition', [100, 100, 1200, 310]);
figure(2)
single_taxis=t_axis(1:segment_duration*sampl_f);
subplot(1,2,1)
plot(single_taxis+1.5,d_sequence(4,:),'r','Linewidth',0.2); hold on;
plot(single_taxis+1.5+2*segment_duration,d_sequence(5,:),'b','Linewidth',0.2)
plot(single_taxis+1.5+segment_duration,d_sequence(idle_index,:),'g','Linewidth',1.5)
xlim([(1/2)*segment_duration+1.5,(5/2)*segment_duration+1.5])
xlabel('Time (s)')
ylabel('Amplitude')
legend('Digit 7','Digit 5','Idle','Location','north');
title('Discrete Time Sequence for keys 7 and 5');

subplot(1,2,2)
plot(single_taxis,d_sequence(4,:),'r','Linewidth',1.2); hold on;
plot(single_taxis,d_sequence(5,:),'b','Linewidth',1.2)
plot(single_taxis,d_sequence(idle_index,:),'g','Linewidth',1.5)
xlabel('Time (s)'); xlim([0,0.012])
ylabel('Amplitude');
title('Discrete Time Sequence for keys 7 and 5');

% savefig(figure(2),'figures/fig3_15.fig')
% saveas(figure(2),'figures/forlatex/fig3_15','epsc')

%% -----------------------------------------------------------------------%
% 3.4.2 - spectrogram


figure (3)
subplot(1,2,1)
spectrogram(y,hanning(Lseg),0,Lseg,sampl_f,'yaxis');
title({'Spectrogram of Landline number '; num2str(line_num)});
ylim([0 2.5])




[S,F] = spectrogram(y,hanning(Lseg),0,Lseg,sampl_f,'yaxis');

subplot(1,2,2)
digits=[1, 2, 4, 5, 6, 7, 11; 0, 2, 7, 5, 8, 3, 6];
plot([F,F,F,F,F,F,F]/1000,20*log10(abs(S(:,2*digits(1,:)-1))),'Linewidth',1.5)
% plot([F,F,F,F,F,F,F]/1000,(abs(S(:,2*digits(1,:)-1))),'Linewidth',1.5)
xlim([0 2.5])
legend(num2str(digits(2,:)'),'location','northeast')
xlabel('Frequency (kHz)')
ylabel('Amplitude (dB)')
title("FFT of Landline' s digit sequences ")

% savefig(figure(3),'figures/fig3_16.fig')
% saveas(figure(3),'figures/forlatex/fig3_16','epsc')

%% -----------------------------------------------------------------------%
% 3.4.4 -Repeat with added noise of low 0.01 medium 1 and high 100 variance

vars=[0.01 0.5 25]';

noise=sqrt(vars).*randn(1,length(y));%generate noise

y_noisy=[y;y;y]+noise;%add the noise


for subplot_index=1:length(vars)
    figure(3+subplot_index)
        
        y=y_noisy(subplot_index,:);
    
        subplot(1,3,1)
        plot(t_axis,y,'Color',RGB1,'Linewidth',0.2);
        title('Landline dual-tone signal');
        xlabel('Time (s)'); xlim([0,21*0.25])
        ylabel('Amplitude')

%         sound(y,sampl_f)


        subplot(1,3,3)
        spectrogram(y,hanning(Lseg),0,Lseg,sampl_f,'yaxis');
        title("Landline' s Spectrogram");
        ylim([0 2.5])

        [S,F] = spectrogram(y,hanning(Lseg),0,Lseg,sampl_f,'yaxis');

        subplot(1,3,2)
        digits=[1, 2, 4, 5, 6, 7, 11; 0, 2, 7, 5, 8, 3, 6];
        plot([F,F,F,F,F,F,F]/1000,20*log10(abs(S(:,2*digits(1,:)-1))),'Linewidth',1.5)
        % plot([F,F,F,F,F,F,F]/1000,(abs(S(:,2*digits(1,:)-1))),'Linewidth',1.5)
        xlim([0 2.5])
        legend(num2str(digits(2,:)'),'location','northeast')
        xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB)')
        title("FFT of digit sequences ")
        


end

%         savefig(figure(4),'figures/fig3_17.fig')
%         saveas(figure(4),'figures/forlatex/fig3_17','epsc')
%         
%         savefig(figure(5),'figures/fig3_18.fig')
%         saveas(figure(5),'figures/forlatex/fig3_18','epsc')
%         
%         savefig(figure(6),'figures/fig3_19.fig')
%         saveas(figure(6),'figures/forlatex/fig3_19','epsc')
