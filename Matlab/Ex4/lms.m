%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 - The least mean square(LMS)algorithm  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:
% x : input sequence (dimensions: Nx1)
% z : output with additive noise(eta (h)) (assumed dimensions: Nx1)
% mu : adaptation gain
% f_order : filter's order (=Nw+1)
% Output:
% y_hat : output of adaptive filter
% e : error
% w_evolution : adaptive weight evolution  (Nw+1)xN matrix

function [y_hat, e, w_evolution] = lms(x, z, mu, f_order)

    N=length(x);
    
    w=zeros(f_order,1); %weights array
    w_evolution = zeros(f_order,N);
    
    for n = f_order:N %since for n<0 x=0
        
        x_segment=x(n:-1: n-f_order+1);
        y_hat(n)=w' *x_segment;
        
        e(n)=z(n)-y_hat(n);
        
        w = w+mu*x_segment*e(n);
        w_evolution(:,n)=w;
    end

end