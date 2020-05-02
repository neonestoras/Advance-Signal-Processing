%%%%%%%%%%%%%%%%%%%%%%%
% 4.3 - Gear shifting %
%%%%%%%%%%%%%%%%%%%%%%%

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

%Initialising Random number generation seed
% seed = rng; 

%Generate a 1000-samples WGN
N=1000;
% x=randn(N,1);
load 'part4_1000wgn.mat';
varx=var(x);


%unknown system filter coefficients
b=[1,2,3,2,1];
a=1;

%output of unknown filter with wgn input x
y_out=filter(b,a,x);
stdy_out=std(y_out);
y=detrend(y_out)/stdy_out;


%Generating WGN h (added noise)
h = (0.1).*randn(1000,1);
stdh = std(h);
varh = var(h);


%realistic output z
z=y+h;
varz=var(z);


%% -----------------------------------------------------------------------%

Nw=4;
f_order=Nw+1;
mu=0.01;
%original lms
[y_hat, e, w_evo]=lms(x,z,mu,f_order);
w_evo=w_evo*stdy_out;

% Gear shifting LMS
mu_upperbound=1/(varx*f_order); %Mean convergence critical mu value
overshoot_lim= (1.2)*b';
[y_hat_gs, e_gs, w_evo_gs, mu_evo] = gs_lms(x, z, mu, mu_upperbound, f_order, overshoot_lim./stdy_out );
w_evo_gs=w_evo_gs*stdy_out;


Lx=length(x);
n=[0:Lx-1]; %discrete time axis

figure (1)
subplot(1,3,1)
plot(n,e.^2,'g','Linewidth',1.3) ; hold on;
plot(n,e_gs.^2,'r','Linewidth',1.1) ;
title('Evolution of squared error') ;
xlabel('Discrete time')
ylabel('Squared error value')
legend('Original LMS','Gear shifting LMS','Location','northeast','fontsize',12)
set(gca,'Color','k')
set(legend, 'Color', 'w')


subplot(1,3,2)
plot(n,w_evo,'Linewidth',1.5);
ylim([0,4]);
xlabel('Discrete time')
ylabel('Magnitude')
title({'Coefficient Evolution', ['for constant \mu= ',num2str(mu)]})

subplot(1,3,3)
plot(n,w_evo_gs,'Linewidth',1.5);hold on;
for i=1:f_order
    yline(overshoot_lim(i),'k--');
end
ylim([0,4]);
xlabel('Discrete time');
ylabel('Magnitude')
title({'Coefficient Evolution', 'with gear shifting'})
legend('w_{opt}[1]','w_{opt}[2]','w_{opt}[3]','w_{opt}[4]','w_{opt}[5]','max(w)');
legend('Fontsize',10)






%% measure the rise time
for i=1:f_order
    or_10(i)=min(n(w_evo(i,:)>0.1*b(i)));
    gs_10(i)=min(n(w_evo_gs(i,:)>0.1*b(i)));
    or_90(i)=min(n(w_evo(i,:)>0.9*b(i)));
    gs_90(i)=min(n(w_evo_gs(i,:)>0.9*b(i)));
end

risetime(:,1)=or_90-or_10;
risetime(:,2)=gs_90-gs_10;

mean_rT=mean(risetime,1);

Type=[{'original'; 'Gear Shifting'}];
Risetime=mean_rT';
recorded_risitimes=table(Type,Risetime)

