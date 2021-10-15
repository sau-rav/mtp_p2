% Gradient Descent method to solve 
% solve min_{T} || AT - TL || such that T PD 
% 
% ****** Input ******
% A, B_inv, B, T, L : five matrices
% S  : initialialization for variable S 
% maxiter: max number of iterations (default = 1e4)
% ****** Output ******
% T   : T PD 

function [T, e, t] = gradDescS(A, B, B_inv, T, L)
    n = size(A, 1); 
    maxiter = 1000; 
    cput = cputime; 
    timemax = Inf;

    e(1) = 0.5 * norm((eye(n) - B * B_inv) * (A * T - T * L) , 'fro')^2; 
    t(1) = cputime - cput;  

    % Step Length
    AAt = A' * A;
    LLt = L * L';
    c = eye(n) - B * B_inv;
    LC = eigs(c, 1) * (eigs(AAt, 1) + eigs(LLt, 1) + 2 * eigs(L, 1) * eigs(A, 1));

    % parameter for the line search
    lsparam = 2; 
    inneriter = 20;
    step = 1 / LC; 
    
    i = 1;  
    while i <= maxiter && cputime-cput <= timemax
        gT = norm(eye(n) - B * B_inv, 'fro')^2 * (A' * (A * T - T * L) - (A * T - T * L) * L'); 
        
        e(i+1) = +Inf; 
        inner = 0; 
        while e(i+1) > e(i) && inner <= inneriter
            Tn = T - gT * step; 
            Tn = projectPSD(Tn); 
            
            e(i+1) = 0.5 * norm((eye(n) - B * B_inv) * (A * Tn - Tn * L) , 'fro')^2; 
            t(i+1) = cputime-cput; 
            
            step = step / lsparam;
            inner = inner + 1; 
        end
        if inner >= inneriter   
           return; 
        end
        T = Tn;
        i = i + 1;
        step = 1 / LC; 
    end
end