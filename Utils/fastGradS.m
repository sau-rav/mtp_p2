% Fast Gradient method to solve 
% solve min_{T} || AT - TL || such that T PD 
% 
% ****** Input ******
% A, B_inv, B, T, L : five matrices
% S  : initialialization for variable S 
% maxiter: max number of iterations (default = 1e4)
% ****** Output ******
% T   : T PD 

function T = fastGradS(A, B, B_inv, T, L) 
    n = size(A, 1); 
    % step size
    AAt = A' * A;
    LLt = L * L';
    c = eye(n) - B * B_inv;
    LC = eigs(c, 1) * (eigs(AAt, 1) + eigs(LLt, 1) + 2 * eigs(L, 1) * eigs(A, 1));
    
    alpha0 = 0.1; 
    alpha(1) = alpha0;
    Yt = T;
    
    i = 1; delta = 1e-3; 
    eps = 1; eps0 = 0; 
    maxiter = 1000;

    while i <= maxiter && eps >= delta * eps0
        Tp = T; 
        alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
        beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
         
        grad = norm(eye(n) - B * B_inv, 'fro')^2 * (A' * (A * T - T * L) - (A * T - T * L) * L'); 
        
        T = projectPSD(Yt - grad / LC);
        Yt = T + beta(i) * (T - Tp); 
        
        if i == 1
            eps0 = norm(Tp - T, 'fro') / (norm(T, 'fro') + 1e-6); 
        end
        eps = norm(Tp - T, 'fro') / (norm(T, 'fro') + 1e-6); 
        i = i + 1; 
    end  
end