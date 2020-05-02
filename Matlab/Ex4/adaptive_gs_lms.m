%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAPTIVE LMS algorithm with Gear Shifting     %
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

function [x_hat, e, w_evolution, mu_evolution] = adaptive_gs_lms(x, mu_begin, mu_max, f_order)
    capped_counter=0;

    N=length(x);
    mu_minimum=0.002;
    decay_constant=2.1;
    
    gearshift_coeff=[mu_max (mu_minimum-mu_max) -1/decay_constant];
    
    w=zeros(f_order,1); %weights array
    w_evolution = zeros(f_order,N);
    mu_evolution= zeros(1,N);
    
    
    mu_evolution(f_order+1)=mu_begin;
    
    
    for n = f_order+1:N %since for n<0 x=0
        
        x_segment=x(n-1:-1: n-f_order);
        x_hat(n)=w' *x_segment;
        
        e(n)=x(n)-x_hat(n); %calculate error
        
        w = w+mu_evolution(n)*x_segment*e(n);
        
        
        w_evolution(:,n)=w;
        
        % ----Gear shift----
        
        % calculate the next mu
        mu_next=gearshift_coeff(1:2)*[1; exp(gearshift_coeff(3)*e(n)^2)];
        
        mu_evolution(n+1)=mu_next;
        
    end
end