
%%  1.1 Statistical estimation  %%
%--------------------------------%
clc;
clear all;
close all;

%Set default sizes
set(groot, 'defaultFigurePosition', [100, 100, 1200, 300]);
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLegendFontSizeMode', 'manual');
%Remove extra whitespace around figures
set(groot,'defaultAxesLooseInset',[0,0,0,0]);
%Show grid on figures
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');


signals=zeros(1000,2);
signals(:,1)=rand(1000,1);%for up to part 4
%==============================================
%==============================================
signals(:,2)=randn(1000,1);%part5 repetition
for p=1:2
%==============================================
%==============================================
    
    
    %%----------------------------------------------------
    %create a signal x[n] of 1000 realisations x=rand(1000,1)
    %that where each x[n] is a realisation of a 
    %uniform random variable X ? U(0, 1) at time instant n

    x=signals(:,p);

    figure(1+(p-1)*4)
    plot(x,'y');set(gca,'Color','k');
    title('Realisation of x[n] for 1000 samples');
    xlabel('Sample number');ylabel('Amplitude');

    %% Part1------------------------------------------------------
    %compute the sample mean
    sample_mean = mean(x)
    mean_bias=sample_mean-0.5
    %calculated mean of the uniform
    %distribution from 0 to 1 is 0.5
    %mean percentage error
    mean_error_percentage=abs(mean_bias/0.5)*100

    %% Part2----------------------------------------------------------
    %Computing the sample standard deviation of x
    sample_std=std(x)
    sample_var_uni=var(x)
    %std percentage error calculated 
    %variance of the uniform distribution [0,1] is 1/12
    std_error=abs(sample_std-sqrt(1/12))
    std_error_percentage=100*std_error/(sqrt(1/12)) %error as a percentage
    %% Part3----------------------------------------------------------------
    %Generate an enseble of 10 1000-sample realisations of X
    %then calculate  the sample means and standard deviations for each realisation
    %and plot

    ensemble_x =rand(1000,10);
    %calculate the mean and std of each column
    x_ens_means=mean(ensemble_x);
    x_ens_std=std(ensemble_x);

    figure(2+(p-1)*4)
    %plotting the means around the theoretical value of 0.5(red line)
    subplot(1,2,1);stem(x_ens_means,'Linewidth',1,'Color','b');hold on
    %stem graph is chosen because the x axis is descrete independent values and
    %the should not be connected together 
    line([1,10],[0.5,0.5],'Color','red');%theoretical value line
    title('Plot of sample means for 10 realisations');
    xlabel('Realisation Index'); ylabel('Sample Mean');
    xlim([1,10]);
    legend('Sample Mean','Theoretical mean','Location','southeast');
    hold off;


    %plotting the stds around the theoretical value
    subplot(1,2,2); stem(x_ens_std,'Linewidth',1,'Color','b'); hold on;
    line([1,10],[1/sqrt(12),1/sqrt(12)],'Color','red','Linewidth',1);%theoretical value indication
    xlabel('Realisation Index'), ylabel('Sample St.Deviation');
    title('Plot of sample standard deviations for 10 realisations');
    legend('Sample STD','Theoretical STD','Location','southeast');
    xlim([1,10]);ylim([0,0.4]);
    hold off

    %% Part4----------------------------------------------------------------
    % Approximate the pdf of X by showing in the same plot the histogram of x,
    % normalised by the number of samples considered,
    % and the theoretical pdf.
    num_ofsamples=[100, 250, 500, 1000, 10000];
    num_ofbins=[10, 25, 50, 100, 1000];
    %effect 1 - of number of generated samples

    figure(3+(p-1)*4)
    for subplot_index=1:length(num_ofsamples)

       subplot(1,length(num_ofsamples),subplot_index);
       if p==1 
           line1=line([0,1],[1,1],'Color','red','Linewidth',1.2);
           hold on;
           xlim([0,1]);ylim([0,1.5]);
           x=rand(num_ofsamples(subplot_index),1);
       elseif p==2
           norm_distr_samples=-4:0.001:4;
           line1=plot(norm_distr_samples,normpdf(norm_distr_samples,0,1),'r','Linewidth',1.2);
           hold on;
           xlim([-4,4]);
           x=randn(num_ofsamples(subplot_index),1);
       end
       histogram(x,num_ofbins(1),'Normalization','pdf');
       title([num2str(num_ofsamples(subplot_index)) ' samples']);
       xlabel('x'); ylabel('pdf(x)');
       uistack(line1,'top');
       hold off;
    end
   
    %effect 2 - of histogram bin numbers
    figure(4+(p-1)*4)
    for subplot_index=1:length(num_ofbins)

       subplot(1,length(num_ofbins),subplot_index);
       if p==1
           line1=line([0,1],[1,1],'Color','red','Linewidth',1.2);
           hold on;
           xlim([0,1]); ylim([0,2.3]);
           x=rand(num_ofsamples(5),1);
       elseif p==2
           norm_distr_samples=-4:0.001:4;
           line1=plot(norm_distr_samples,normpdf(norm_distr_samples,0,1),'r','Linewidth',1.2);
           hold on;
           xlim([-4,4]);
           x=randn(num_ofsamples(5),1);
       end
       histogram(x,num_ofbins(subplot_index),'Normalization','pdf');
       title([num2str(num_ofbins(subplot_index)) ' bins']);
       xlabel('x'); ylabel('pdf(x)');
       uistack(line1,'top');
       hold off;
    end

%% Part5----------------------------------------------------------------
end %end repeting parts 1 to 4 with x=randn(1000,1)


%% change of mean and std as function of sample size
N=10000;
sample_m=zeros(N,2)
sample_stds=zeros(N,2)
for s=1:1:N
    unif=rand(s,1);
    gauss=randn(s,1);
    sample_m(s,1)=mean(unif);
    sample_m(s,2)=mean(gauss);
    sample_stds(s,1)=std(unif);
    sample_stds(s,2)=std(gauss);
end

figure(9)
subplot(1,2,1)%sample mean
samples_axis=[1:N];
%uniform distribution mean
plot(samples_axis,sample_m(:,1),'Linewidth',1,'Color','g');hold on;
line([1,N],[0.5,0.5],'Color','red','Linewidth',1.5);
set(gca,'Color','k');
%gaussian distribution mean
plot(samples_axis,sample_m(:,2),'y','Linewidth',1);
line([1,N],[0,0],'Color','b','Linewidth',1.5);
xlabel('\fontsize{14}sample size'), ylabel('\fontsize{14}Magnitude of sample mean')
ylim([-0.2 0.6]);
title('\fontsize{14}\bfMean vs sample size')
legend('\fontsize{10}sample mean for uniform distr.','\fontsize{10}theoretical mean for uniform distr.','\fontsize{10}sample mean for gaussian distr.','\fontsize{10}theoretical mean for gaussian distr.','Location','east');
set(legend,'color','w');
hold off;

subplot(1,2,2)
%uniform distribution std
plot(samples_axis,sample_stds(:,1),'g','Linewidth',1);hold on;
line([1,N],[(sqrt(1/12)),(sqrt(1/12))],'Color','red','Linewidth',1.3);
%gaussian distribution std
plot(samples_axis,sample_stds(:,2),'y','Linewidth',1);
line([1,N],[1,1],'color','b','Linewidth',1.5);

set(gca,'Color','k');
xlabel('\fontsize{14}sample size'), ylabel('\fontsize{14}Magnitude of sample std')
ylim([0.2 1.2]);
title('\fontsize{14}\bfStandard deviation vs sample size')
legend('\fontsize{10}sample std for uniform distr.','\fontsize{10}theoretical std for uniform distr.','\fontsize{10}sample std for gaussian distr.','\fontsize{10}theoretical std for gaussian distr.','Location','east');
set(legend,'color','w');
hold off;