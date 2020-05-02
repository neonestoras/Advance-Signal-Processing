%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.3 The Least Squares Estimation(LSE)of AR Coeffcients  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 350]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);


%% -----------------------------------------------------------------------%
% 3.3.3 find the coefficients a using LSE approach
load sunspot.dat; %Data obtained from MATLAB itself
clc
suns=zeros(288,2);
suns(:,1)=sunspot(:,2); %original data
suns(:,2)=(sunspot(:,2)-mean(sunspot(:,2)))./std(sunspot(:,2)); %standardised data
x=suns(:,2);


M=length(x);
av_coefficient_diff=zeros(1,10);
e_yw=zeros(1,10);
lse_error=zeros(1,10);
figure(1)
for P=1:10
   
    %create the observation matrix H
    H=zeros(M,P);
    for m=1:1:M
        for p=1:1:P
            if m-p>0
                H(m,p) = x(m-p);
            else
                break
            end
        end
        
    end
    
    % LSE approach (estimation)
    a_lse=(inv((H'*H)))*H'*x;
    lse_error(P)=((x-H*a_lse)'*(x-H*a_lse))/M;
    corrected_lse_error(P)=lse_error(P)+(((2*P)*(P+1))/(M-P-1));
    
    [h_lse(:,P),w_lse(:,P)]= freqz(1,[1 -a_lse'],M); %create the AR model
    
    
    %Yule-Walker approach
    [a_yw,e_yw(P)] = aryule(x,P);
    a_yw=-a_yw(2:length(a_yw)); %remove the first element that is always 1 and negative
    corrected_ym_error(P)=e_yw(P)+(((2*P)*(P+1))/(M-P-1));
    
    
    %difference in coeffiients
    av_coefficient_diff(P)=mean(abs(a_lse-a_yw'));
    
    
    
    
    %plot
    subplot(2,5,P)
    power_index=1:1:P;
    stem(power_index,a_yw,'b','Linewidth',1.6); hold on;
    stem(a_lse,'r', 'MarkerFaceColor','red','MarkerEdgeColor','none','Linewidth',1.2);
    title(['Model Order ',num2str(P)]);
    xlabel('Coefficient index');
    xlim([0,P+1]); xticks(1:P);
    ylabel('Amplitude');
    
end
legend('Yule-Walker','LSE')
legend('Location', 'southoutside', 'Fontsize',11)


%% 3.3.4- Approximation error
%%-----------------------------------------------------------------------%
set(groot, 'defaultFigurePosition', [100, 100, 1100, 300]);
figure(2)
p=1:1:P;
subplot(1,2,2)
plot(p,corrected_lse_error,'r','Linewidth',2);hold on;
plot(p,corrected_ym_error,'b--','Linewidth',1.8);
title('Corrected Error')
ylabel('Magnitude (AU)')
xlabel('Model Order')
legend('\fontsize{10}LSE','\fontsize{10}Yule-Walker','Location','north')
xlim([1,10])

subplot(1,2,1)
plot(p,(lse_error),'r','Linewidth',2);hold on;
plot(p,(e_yw),'b--','Linewidth',1.8);
title('Estimation Error')
ylabel('Magnitude (AU)')
xlabel('Model Order')
legend('\fontsize{10}LSE','\fontsize{10}Yule-Walker','Location','north')
xlim([1,10])

% subplot(1,2,1)
% stem(p,(av_coefficient_diff),'Linewidth',2, 'Color',[0.4660 0.6740 0.1880]);
% title('Coefficient difference')
% ylabel('Difference (AU)')
% xlabel('Model Order')
% xlim([1,10])

%% 3.3.5 - PSD of the obtained models using LSE
%%-----------------------------------------------------------------------%
[pgm_x freq_x]=pgm(x');
set(groot, 'defaultFigurePosition', [100, 100, 1200, 350]);
figure(3)
for p=1:1:P
    subplot(2,5,p)
        plot(freq_x,pgm_x,'k','Linewidth',0.8);hold on;
    
        plot(w_lse(:,p)/(2*pi),lse_error(p)*(abs(h_lse(:,p)).^2),'b','Linewidth',1.25);
        xlabel('\fontsize{10}Normalised Frequency');
        ylabel('Magnitude');
        str=sprintf('LSE, AR(%d) model',p);
        title(str)
        xlim([0 0.5]);
end
legend('\fontsize{9}Periodogram','\fontsize{9}LSE-Based model','Location','northwest');

%% 3.3.6 -  behaviour of MSE versus N for model optimal order (2)
% and find optimum data length N
%%-----------------------------------------------------------------------%

P=2;
N=10:5:250;
for subplot_index=1:length(N)
   n=N(subplot_index);
   seq=x(1:n);
   
    %create the observation matrix H
    H=zeros(n,P);
    for m=1:1:n
        for p=1:1:P
            if m-p>0
                H(m,p) = seq(m-p);
            else
                break
            end
        end
        
    end
    
    % LSE approach (estimation)
    a_lse=(inv((H'*H)))*H'*seq;
    lse_error(subplot_index)=((seq-H*a_lse)'*(seq-H*a_lse))/n;
    
end

set(groot, 'defaultFigurePosition', [100, 100, 500, 260]);
figure (4)
plot(N, lse_error,'g','Linewidth',2)
ylabel('MSE')
xlabel('Data length')
title('LSE model MSE with data-length')
set(gca,'Color','k')

% savefig(figure(4),'figures/fig3_13.fig')



