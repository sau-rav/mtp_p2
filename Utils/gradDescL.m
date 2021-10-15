% Gradient Descent method to solve 
% solve min_{L} || AT - TL || such that ||L|| <= 1 
% 
% ****** Input ******
% A, B_inv, B, T, L : five matrices
% L  : initialialization for variable L 
% maxiter: max number of iterations (default = 1e4)
% ****** Output ******
% L   : ||L|| <= 1 

function [L, e, t] = gradDescL(A, B, B_inv, T, L)   
    n = size(A, 1); 
    maxiter = 1000; 
    cput = cputime; 
    timemax = Inf;
    
    e(1) = 0.5 * norm((eye(n) - B * B_inv) * (A * T - T * L) , 'fro')^2; 
    t(1) = cputime - cput;  

    % Step Length  
    TTt = T' * T;
    c = eye(n) - B * B_inv;
    LC = eigs(c, 1) * eigs(TTt, 1);

    % parameter for the line search
    lsparam = 2; 
    inneriter = 20;
    step = 1 / LC; 
    
    i = 1;  
    while i <= maxiter && cputime-cput <= timemax
        gL = -norm(eye(n) - B * B_inv, 'fro')^2 * T' * (A * T - T * L); 
        
        e(i+1) = +Inf; 
        inner = 0; 
        while e(i+1) > e(i) && inner < inneriter 
            
            Ln = L - gL * step; 
            Ln = projectNorm(Ln); 
            
            e(i+1) = 0.5 * norm((eye(n) - B * B_inv) * (A * T - T * Ln) , 'fro')^2; 
            t(i+1) = cputime-cput; 
            
            step = step / lsparam;
            inner = inner + 1; 
        end
        if inner >= inneriter   
           return; 
        end
        L = Ln;
        i = i + 1; 
        step = 1 / LC; 
    end
end