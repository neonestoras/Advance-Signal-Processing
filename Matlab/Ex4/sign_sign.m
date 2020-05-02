%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sign Error Algorithm  %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:
% x : input sequence (dimensions: Nx1)
% mu : adaptation gain
% f_order : filter's order (=Nw+1)
% Output:
% y_hat : output of adaptive filter
% e : error
% w_evolution : adaptive weight evolution  (Nw+1)xN matrix

function [x_hat, e, w_evolution] = sign_sign(x, mu, f_order)

    N=length(x);
    
    w=zeros(f_order,1); %weights array
    w_evolution = zeros(f_order,N);
    
    x_hat=[];
    e=[];
    
    for n = f_order+1:N %since for n<0 x=0
        
        x_segment=x(n-1:-1: n-f_order);
        x_hat(n)=w' *x_segment;
        
        e(n)=x(n)-x_hat(n);
        
        w = w+ mu * sign(x_segment) *sign(e(n));
        w_evolution(:,n)=w;
    end

end