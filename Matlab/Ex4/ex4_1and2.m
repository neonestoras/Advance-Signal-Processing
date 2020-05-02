%%%%%%%%%%%%%%%%%%%%%%%
% 4.1 - Wiener filter %
%%%%%%%%%%%%%%%%%%%%%%%

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


%% -----------------------------------------------------------------------%
% 4.1.1 - Rxx and pxx statistics

%Generate a 1000-samples WGN
N=1000;
% x=randn(N,1);
load 'wgn1000.mat'


%unknown system filter coefficients
b=[1,2,3,2,1];
a=1;

%output of unknown filter with wgn input x
y_out=filter(b,a,x);
stdy_out=std(y_out);
y=detrend(y_out)/stdy_out;


%Generating WGN h (added noise)
% h = (0.1).*randn(1000,1);
load '0p1_eta.mat'
stdh = std(h);
varh = var(h);


%realistic output z
z=y+h;
varz=var(z);
snrz=snr(y,h);% Calculate the Signal-toNoise Ratio (SNR) in dB for z[n].


% the ACF
Nw = 4; % max lag used
[acfx, lags_xx]=xcorr(x,Nw,'unbiased');

[acf_xz, lags_xz] = xcorr(x,z,Nw,'unbiased');
pzx=flip(acf_xz(lags_xz<=0));

Rxx=zeros(Nw+1);
for column=0:-1:-Nw
   Rxx(:,-column+1)=acfx(lags_xx>=column & lags_xx<=Nw+column);
end

%optimal coefficients matrix
wopt=(inv(Rxx)*pzx)*stdy_out;


%% -----------------------------------------------------------------------%
% 4.1.2 - testing different variance additive noise (h)

vars=[0.1,0.5,1,3,6,10];
stds=sqrt(vars);


%generate additive noice with different variance
load 'variance_test_etas.mat'
% h2=zeros(1000,length(vars));
% for i=1:length(vars)
%     h2(:,i)=stds(i)*randn(1000,1);
% end

y2=y*ones(1,length(vars));

z2=h2+y2;


snr_z2=zeros(1,6);
for i=1:length(vars)
    snr_z2(i)=snr(y2(:,i),h2(:,i));
    [corr_xz, lags_xz] = xcorr(x,z2(:,i),Nw,'unbiased');
    pzx=flip(corr_xz(lags_xz<=0));
    wopt2(:,i)=(inv(Rxx)*pzx)*stds(i);
end

theoretical_snr=round((10*log10(vars))',4);
measured_snr=round((snr_z2)',4);
obtained_wopt=round((wopt2)',4);

VarianceEffect=table(vars',theoretical_snr, measured_snr,obtained_wopt );

%% Repeat for Nw > 4
Nw_tested=[4:7]';
wopt3=zeros(Nw_tested(length(Nw_tested))+1,length(Nw_tested)-1);
for j=2:length(Nw_tested)
    
    Nw=Nw_tested(j);
    
    [acfx, lags_xx]=xcorr(x,Nw,'unbiased');

    [acf_xz, lags_xz] = xcorr(x,z,Nw,'unbiased');
    pzx=flip(acf_xz(lags_xz<=0));

    Rxx=zeros(Nw+1);
    for column=0:-1:-Nw
       Rxx(:,-column+1)=acfx(lags_xx>=column & lags_xx<=Nw+column);
    end

    %optimal coefficients matrix
    w_temp=(inv(Rxx)*pzx)*stdy_out;
    
    %make all w vsctors of the same length
    for i=Nw+1:Nw_tested(length(Nw_tested))
        w_temp=[w_temp;0];
    end
   
    wopt3(:,j-1)=w_temp;

end

wopt_obtained=round(([[wopt;0;0;0],wopt3])',4);

coeff_number_effet=table(Nw_tested,wopt_obtained);




%% -----------------------------------------------------------------------%

%% -----------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 - The least mean square(LMS)algorithm  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -----------------------------------------------------------------------%
% 4.2.1 - test the lms function 

Nw=4;
f_order=Nw+1;


% Calculate the critical mu value
E = eig(Rxx);
mean_convergence_crit_mu=2/max(E)
ms_convergence_crit_mu=2/trace(Rxx)

mu=0.003;

[y_hat, e, w_evo]=lms(x,z,mu,f_order);
w_evo=w_evo*stdy_out;

Lx=length(x);
n=[0:Lx-1]; %discrete time axis
zero_axis=zeros(Lx,1);

figure(1)
subplot(1,2,1)
plot(n,z,'y','Linewidth',1.2);hold on;
plot(n,y_hat,'b','Linewidth',1.3);
% stem3(n,zero_axis,e,'.')
xlabel('\fontsize{13}Discrete Time (n)')
ylabel('\fontsize{13}Amplitude (AU)')
set(gca,'Color','k')
% legend('True Output', 'Estimated Output', 'Estimation error')
legend('True Output', 'Estimated Output')
legend('Location','south','orientation','horizontal','Fontsize',12)
set(legend,'color','w');

subplot(1,2,2)
plot3([Lx;Lx;Lx;Lx;Lx],zeros(length(b),1),b','go','Linewidth',5);hold on;
plot3([n;n;n;n;n],zero_axis,w_evo,'Linewidth',1.5);
plot3(n,e,zero_axis);
ylabel('\fontsize{13}Amplitude(AU)')
zlabel('\fontsize{13}Amplitude(AU)')
xlabel('\fontsize{13}Discrete Time (n)')
legend('\fontsize{10}Actual weights','\fontsize{10}Estimated','\fontsize{10}weights','.','.','.','\fontsize{10}Estimation Error')


%% -----------------------------------------------------------------------%
% 4.2.2 - adaptation gain mu=0.01
mu=0.01;

[y_hat, e, w_evo]=lms(x,z,mu,f_order);
w_evo=w_evo*stdy_out;



figure (2)
subplot(1,2,1)
plot(n,e.^2,'g','Linewidth',1.3) ;
title('Evolution of squared error over time') ;
xlabel('Discrete time')
ylabel('Squared error value')
set(gca,'Color','k')


subplot(1,2,2)
plot(n,w_evo,'Linewidth',1.5);
xlabel('Discrete time')
ylabel('Magnitude')
title('Evolution of Filter coefficients over time')
legend('w_{opt}[1]','w_{opt}[2]','w_{opt}[3]','w_{opt}[4]','w_{opt}[5]');
legend('Fontsize',11)


%% effect of adaptive gain

adaptive_gains=[0.002 0.02 0.2 0.4];

figure (3)

for m=1:length(adaptive_gains)
    
    subplot(1,length(adaptive_gains),m)
    
    mu=adaptive_gains(m);
    
    [y_hat, e, w_evo]=lms(x,z,mu,f_order);
    w_evo=w_evo*stdy_out;
    
    plot(n,w_evo,'Linewidth',1.5);
    xlabel('Discrete time')
    ylabel('Magnitude')
    title({'Coefficient Evolution', ['for \mu = ',num2str(mu)]})
    
end
legend('w_{opt}[1]','w_{opt}[2]','w_{opt}[3]','w_{opt}[4]','w_{opt}[5]');
legend('Fontsize',11)





