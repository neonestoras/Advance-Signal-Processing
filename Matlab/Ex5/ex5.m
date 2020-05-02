%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MLE for the Frequency of a signal%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
N=10;
n=0:1:N-1;
tested_f=[0.05 0.2 0.3 0.45];
figure (1)
for f_i=1:4
    f=tested_f(f_i);
    x=cos(2*pi*f*n);
    
%     [pgm_x, freq]=pgm(x);
    [pgm_x, freq]=periodogram(x);
    [~,I_max]=max(pgm_x);
    
    mle(f_i)=freq(I_max)/(2*pi);
    
    subplot(1,4,f_i)
    plot(freq./(2*pi),pgm_x,'b');hold on;
    xline(mle(f_i),'r--');
    ylabel('Magnitude')
    xlabel('Normalised Frequency')
    title([ 'f_{0} = ',num2str(f)]);
    xlim([0 0.5])
    
end
legend('\fontsize{10}Periodogram','\fontsize{10}MLE')
% savefig(figure(1),'figures/fig5_1.fig')

%%
tested_frequencies=tested_f';
measured_mle=mle';
error=(abs(mle-tested_f)./tested_f)';

mle_table=table(tested_frequencies,measured_mle,error);




