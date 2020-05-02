%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.6 Dealing with computational complexity: signalgorithms %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%Synthesis
    N=1000;
    h=randn(N,1);

    % AR model coefficients
    a=[1 0.9 0.2];
    b=1;

    % Filter h WN
    x=filter(b,a,h);
    
    p=length(a)-1; %order
    mu=0.02; %adaptation gain

%% Adaptive LMS
[y, e,w] = adaptive_lms(x,mu,p);

%% Signed-error
[y_se, e_se,w_se] = signed_error(x,mu,p);

%% Signed-regressor
[y_sre, e_sre,w_sre] = signed_regressor(x,mu,p);

%% Sign - Sign
[y_ss, e_ss,w_ss] = sign_sign(x,mu,p);
 

%% -----------------------------------------------------------------------%
Lx=length(x);
n=[0:Lx-1]; %discrete time axis
% 
% figure (1)
% subplot(2,4,1)
%     polarplot(2*pi*n/Lx,w(1,:)-a(2)*ones(1,Lx),'b','Linewidth',1);
%     title('Adaptive LMS')
%     thetaticklabels(['']);
% subplot(2,4,2)   
%     polarplot(2*pi*n/Lx,w_se(1,:)-a(2)*ones(1,Lx),'b','Linewidth',1);
%     title('Signed-Error')
%     thetaticklabels(['']);
% subplot(2,4,3)
%     polarplot(2*pi*n/Lx,w_sre(1,:)-a(2)*ones(1,Lx),'b','Linewidth',1);
%     title('Signed-Regressior')
%     thetaticklabels(['']);
% subplot(2,4,4)
%     polarplot(2*pi*n/Lx,w_ss(1,:)-a(2)*ones(1,Lx),'b','Linewidth',1);
%     title('Sign-Sign')
%     thetaticklabels(['']);
%     legend('a_1')
% 
% subplot(2,4,5)
%     polarplot(2*pi*n/Lx,w(2,:)-a(3)*ones(1,Lx),'r','Linewidth',1);thetaticklabels(['']);
% subplot(2,4,6)   
%     polarplot(2*pi*n/Lx,w_se(2,:)-a(3)*ones(1,Lx),'r','Linewidth',1);thetaticklabels(['']);
% subplot(2,4,7)
%     polarplot(2*pi*n/Lx,w_sre(2,:)-a(3)*ones(1,Lx),'r','Linewidth',1);thetaticklabels(['']);
% subplot(2,4,8)
%     polarplot(2*pi*n/Lx,w_ss(2,:)-a(3)*ones(1,Lx),'r','Linewidth',1);thetaticklabels(['']);
%     legend('a_2')
% 
% savefig(figure(1),'figures/fig4_9.fig')
% saveas(figure(1),'figures/forlatex/fig4_9','epsc')
    
figure (2)
subplot(1,4,1)
    polarplot(2*pi*n/Lx,(-w(2,:)./a(3)-1).*ones(1,Lx),'r','Linewidth',1);hold on;
    polarplot(2*pi*n/Lx,(-w(1,:)./a(2)-1).*ones(1,Lx),'b','Linewidth',1);
    title('Adaptive LMS')
    thetaticklabels(['']);
subplot(1,4,2)   
    polarplot(2*pi*n/Lx,(-w_se(2,:)./a(3)-1).*ones(1,Lx),'r','Linewidth',1);hold on;
    polarplot(2*pi*n/Lx,(-w_se(1,:)./a(2)-1).*ones(1,Lx),'b','Linewidth',1);
    title('Signed-Error')
    thetaticklabels(['']);
subplot(1,4,3)

    polarplot(2*pi*n/Lx,(-w_sre(2,:)./a(3)-1).*ones(1,Lx),'r','Linewidth',1);hold on;
    polarplot(2*pi*n/Lx,(-w_sre(1,:)./a(2)-1).*ones(1,Lx),'b','Linewidth',1);
    title('Signed-Regressor')
    thetaticklabels(['']);
subplot(1,4,4)
    polarplot(2*pi*n/Lx,(-w_ss(2,:)./a(3)-1).*ones(1,Lx),'r','Linewidth',1);hold on;
    polarplot(2*pi*n/Lx,(-w_ss(1,:)./a(2)-1).*ones(1,Lx),'b','Linewidth',1);
    title('Sign-Sign')
    thetaticklabels(['']);
legend('a_2','a_1')

savefig(figure(2),'figures/fig4_9.fig')

figure (3)
subplot(1,4,1)
    plot([n]',-w(1,:)','b');hold on;
    plot([n]',-w(2,:)','r');
    title('Adaptive LMS')
    xlabel('Discrete time')
    ylabel('Magnitude')
    yline(a(2),'k--')
yline(a(3),'k--')
subplot(1,4,2)   
    plot([n]',-w_se(1,:)','b');hold on;
    plot([n]',-w_se(2,:)','r');
    title('Signed-Error')
    xlabel('Discrete time')
    ylabel('Magnitude')
    yline(a(2),'k--')
yline(a(3),'k--')
    subplot(1,4,3)
    plot([n]',-w_sre(1,:)','b');hold on;
    plot([n]',-w_sre(2,:)','r');
    title('Signed-Regressor')
    xlabel('Discrete time')
    ylabel('Magnitude')
    yline(a(2),'k--')
yline(a(3),'k--')
subplot(1,4,4)
    plot([n]',-w_ss(1,:)','b');hold on;
    plot([n]',-w_ss(2,:)','r');
    title('Sign-Sign')
    xlabel('Discrete time')
ylabel('Magnitude')
    yline(a(2),'k--')
yline(a(3),'k--')
legend('a_1','a_2','\fontsize{10}Theoretical','\fontsize{10}Values')

savefig(figure(3),'figures/fig4_10.fig')