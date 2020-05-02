%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 - LMS algorithm with Gear Shifting     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:
% x : input sequence (dimensions: Nx1)
% z : output with additive noise(eta (h)) (assumed dimensions: Nx1)
% mu_begin : adaptation gain starting value
% mu_max : adaptation gain maximum value bound
% f_order : filter's order (=Nw+1)
% max_overshoot : maximum allowed weight values (Nw+1)x1 matrix

% Output:
% y_hat : output of adaptive filter
% e : error
% w_evolution : adaptive weight evolution  (Nw+1)xN matrix
% mu_evolution : adaptation gain evolution (Nw+1)xN matrix

function [y_hat, e, w_evolution, mu_evolution] = gsV3_lms(x, z, mu_begin, mu_max, f_order, max_overshoot)
    capped_counter=0;

    N=length(x);
    mu_minimum=0.002;
    decay_constant=2.1;
    
    gearshift_coeff=[mu_max (mu_minimum-mu_max) -1/decay_constant];
    
    w=zeros(f_order,1); %weights array
    w_evolution = zeros(f_order,N);
    mu_evolution= zeros(1,N);
    
    
    mu_evolution(f_order)=mu_begin;
    
    
    for n = f_order:N %since for n<0 x=0
        
        x_segment=x(n:-1: n-f_order+1);
        y_hat(n)=w' *x_segment;
        
        e(n)=z(n)-y_hat(n); %calculate error
        
        w = w+mu_evolution(n)*x_segment*e(n);
        
        % ----Overshoot correction ------
        for w_i=1:f_order
            if (w(w_i)>max_overshoot(w_i))
                w(w_i)=max_overshoot(w_i);
                capped_counter=capped_counter+1;
            end
        end
        
        w_evolution(:,n)=w;
        
        % ----Gear shift----
        gearshift_coeff(1)=(0.05*mu_evolution(n)+0.95*mu_max);
        
        
        % calculate the (proposed) next mu
        mu_next=gearshift_coeff(1:2)*[1; exp(gearshift_coeff(3)*e(n)^2)];
        
        mu_evolution(n+1)=mu_next;
        
    end
end