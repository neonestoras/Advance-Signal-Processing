%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.1 ACF of uncorrelated and correlated sequences%
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

N=1000;
x=randn(1,1000); %1000-sample WGN realisation

%%  2.1.1 - unbiased estimate of the ACF for a 1000-sample realisation of WGN %% 
%------------------------------------------------------------------------------%
load('x.mat');

[acf,timelag] = xcorr(x,'unbiased'); 
% 'unbiased' - scales the raw correlation by 1/(M-abs(lags)).
% 
% returns the length 2*N-1 auto-correlation sequence of x. 
% The zeroth lag of the output correlation
% is in the middle of the sequence, at element N.
%thus N is in [-(N-1),N-1]

figure (1)
%plot of ACF 
subplot(1,3,2)
plot(timelag,acf,'g','Linewidth',1);
set(gca,'color','k');
title({'ACF estimate for WGN','using |\tau| <1000'}); 
ylabel('Correlation, $$\hat{R(\tau)}$$','Interpreter','Latex');
xlabel('Time Lag (\tau)');    
axis([-N+1 N-1 min(acf) 1]);

%plot of theoretical ACF
subplot(1,3,1)
theor_acf=zeros(1,2*N-1); theor_acf(N)=1;
plot(timelag,theor_acf,'b','Linewidth',1.3);
set(gca,'color','k');
title('Theoretical ACF of WGN'); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag (\tau)');    
axis([-N+1 N-1 min(acf) 1]);

%Illustrating symmetry of autocorrelation function
subplot(1,3,3)
sym=abs(acf-flip(acf));
stem(timelag,sym,'r','.','Linewidth',1.3);
title({'Illustrating symmetry', 'of ACF'}); 
ylabel('$$| \ \hat{R(\tau)}-\hat{R(-\tau)} \ |$$','Interpreter','Latex');
xlabel('Time Lag (\tau)');
axis([-N+1 N-1 min(sym) max(sym)]);

%%  2.1.2 - zoom onto the region |?| < 50 and compare with the shape of the ACF for large ? %% 
%------------------------------------------------------------------------------%

figure (2)
subplot(1,3,2)
plot(timelag, acf,'b','Linewidth',1.2);
zoom xon; zoom(20); %zoom by 20X on x-axis
ylabel('Correlation $$\hat{R(\tau)}$','Interpreter','Latex','Linewidth',1.2);
xlabel('Time Lag (\tau)');

%Plot for larger time lag values
timelag_upper=500:length(acf)/2;
figure(2)
for subplot_index=[1,3]
subplot(1,3,subplot_index)
%lower part
    if subplot_index==1
        plot(-timelag_upper, acf(N-timelag_upper),'r','Linewidth',1.2);
        xlim([-(N-1) -500]);
 %upper part
    elseif subplot_index==3
        plot(timelag_upper, acf(N+timelag_upper),'r','Linewidth',1.2);
        xlim([500 N-1]);
    end
        
ylabel('Correlation $$\hat{R(\tau)}$','Interpreter','Latex','Linewidth',1.2);
xlabel('Time Lag (\tau)');

end
% saveas(figure(2),'figures/fig2_2.fig')
% saveas(figure(2),'figures/forlatex/fig2_2','epsc')


%%  2.1.3 - empirical timelag bound %% 
%------------------------------------------------------------------------------%

% maximum allowed error is 3*std
% std of biased estimator is sqrt(N)
% std of the unbiased ACF is plotted and compared to unbiased
std_acf=zeros(1,2*N-1);
for i =1:1:(2*N-1)
    std_acf(i) = 1/sqrt(1000 - abs(timelag(i))) ;
end

allowed_error=3/sqrt(N);

figure (3)
plot(timelag,std_acf); hold on;
line([-(N-1),N-1],[1/sqrt(N),1/sqrt(N)],'Color','red');
line([-(N-1),N-1],[allowed_error,allowed_error],'Color','black');%allowed error is 3*std
title('std of ACF estimate'); 
ylabel('std magnitude');
xlabel('Time Lag (\tau)'); 



high_count=0;
low_count=0;
for err_index=1:1:1000
    if abs(acf(1000+err_index))>allowed_error
        if high_count==1
            break;
        else
            high_count=high_count+1;
        end
    end
    if abs(acf(1000-err_index))>allowed_error
        if low_count==1
            break;
        else
            low_count=low_count+1;
        end
    end
end
tau_bound=err_index;

figure(4)
plot(timelag,acf,'y','Linewidth',1);hold on;
set(gca,'color','k');
line([-(N-1),N-1],[allowed_error,allowed_error],'Color','red','Linewidth',1.3);%allowed error is 3*std
line([-(N-1),N-1],[-allowed_error,-allowed_error],'Color','red','Linewidth',1.3);%allowed error is 3*std
line([tau_bound,tau_bound],[min(acf), 1],'Linewidth',1.3,'Color','blue','LineStyle','--')
line([-tau_bound,-tau_bound],[min(acf), 1],'Linewidth',1.3,'Color','blue','LineStyle','--')
axis([-N+1 N-1 min(acf) 1]);
legend('ACF','Maximum','allowed error','Empirical', 'Time Lag bound')
set(legend,'Color','w')
ylabel('Correlation $$\hat{R(\tau)}$','Interpreter','Latex','Linewidth',1.2);
xlabel('Time Lag (\tau)');


%%  2.1.4 -  WGN vector x and ?lter it by a moving average (MA) ?lter 
%         - with unit coef?cients of order 9%% 
%------------------------------------------------------------------------------%
x=randn(1,1000); %Generating 1000-sample WGN
y=filter(ones(9,1),1,x);  % first component is numerator and denominator coefficients

figure(5)
plot(x); hold on; plot(y);


figure(6)
acf_y = xcorr(y,'unbiased');
stem(-20:20,acf_y(-20+1000:20+1000),'b','Linewidth',1.2);
title('Unbiased ACF estimate for WGN, when filtered by MA filter of order 9','Linewidth',1.5); 
ylabel('Correlation, R(\tau)');
xlabel('Time Lag (\tau)');
savefig(figure(6),'figures/fig2_4.fig')
saveas(figure(6),'figures/forlatex/fig2_4','epsc')
%%
%effect of filter order
clc
f_order=[5,10,15];

figure(7)
for subplot_index=1:length(f_order)
    subplot(1,length(f_order),subplot_index);
    
    y_2=filter(ones(f_order(subplot_index),1),1,x);
    acf_y_2=xcorr(y_2, 'unbiased');
    stem([-20:1:20],acf_y_2(-20+1000:1:20+1000),'r','.','Linewidth',1.15);
    str=sprintf('MA of order %d',f_order(subplot_index));title(str); 
    ylabel('Correlation, R(\tau)');
    xlabel('Time Lag (\tau)');
end
% savefig(figure(7),'figures/fig2_5.fig')
% saveas(figure(7),'figures/forlatex/fig2_5','epsc')





