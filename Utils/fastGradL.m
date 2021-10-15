% Fast Gradient method to solve 
% solve min_{L} || AT - TL || such that ||L|| <= 1 
% 
% ****** Input ******
% A, B_inv, B, T, L : five matrices
% L  : initialialization for variable L 
% maxiter: max number of iterations (default = 1e4)
% ****** Output ******
% L   : ||L|| <= 1 

function L = fastGradL(A, B, B_inv, T, L) 
    n = size(A, 1); 
    % step size
    TTt = T' * T;
    c = eye(n) - B * B_inv;
    LC = eigs(c, 1) * eigs(TTt, 1);
    
    alpha0 = 0.1; 
    alpha(1) = alpha0;
    Yl = L;
    
    i = 1; delta = 1e-3; 
    eps = 1; eps0 = 0; 
    maxiter = 1000;

    while i <= maxiter && eps >= delta * eps0
        Lp = L; 
        alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
        beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
         
        grad = -norm(eye(n) - B * B_inv, 'fro')^2 * T' * (A * T - T * L); 
        
        L = projectNorm(Yl - grad / LC);
        Yl = L + beta(i) * (L - Lp); 
        
        if i == 1
            eps0 = norm(Lp - L, 'fro') / (norm(L, 'fro') + 1e-6); 
        end
        eps = norm(Lp - L, 'fro') / (norm(L, 'fro') + 1e-6);  
        i = i + 1; 
    end  
end