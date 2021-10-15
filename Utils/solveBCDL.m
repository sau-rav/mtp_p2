
function [S, L, K_final, flag] = solveBCDL(A, B, options)
    
    n = size(A, 1);
    B_inv = pinv(B);
    i = 1;
    K_final = [];
    flag = false;
    
    % Random Initialisation
    S = randn(n); S = S * S'; 
    L = rand(n); L = L / (2 * norm(L));
    
    min_norm = norm((eye(n) - B * B_inv) * (S * A - L * S), 'fro');
    fprintf('norm value initial : %2.10f\n', min_norm);
    e = [];
    e(i) = min_norm;
    
    disp('Solving BCD for checking feasibility');
    while i < options.itermax ...
        && (i <= 3 || e(i - 1) - e(i - 2) > 1e-6)
        
        fprintf('iteration : %d\n', i);
        % keep S as constant solve for L
        Ln1 = gradDescL(A, B, B_inv, S, L);
        disp('Ln1 on gradDesc'); 
        disp(norm((eye(n) - B * B_inv) * (S * A - Ln1 * S), 'fro'));
        cvx_begin sdp quiet
            variable Ln(n,n) 
            minimize(norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro'))
            subject to
                norm(Ln, 'fro') <= 1;
        cvx_end
        disp('Ln on cvx'); 
        disp(norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro'));
        
        % accept grad 
        Ln = Ln1;
        
        if norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro') < min_norm
            L = Ln;
            min_norm = norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro');
        end
        
        % keep L as constand solve for S
        cvx_begin sdp quiet
            variable Sn(n,n)
            minimize(norm((eye(n) - B * B_inv) * (Sn * A - L * Sn), 'fro'))
            subject to
                Sn == semidefinite(n); 
        cvx_end 
        if norm((eye(n) - B * B_inv) * (Sn * A - L * Sn), 'fro') < min_norm
            S = Sn;
            min_norm = norm((eye(n) - B * B_inv) * (Sn * A - L * Sn), 'fro');
        end
        
        fprintf('norm value at iteration end : %2.10f\n', min_norm);
        if norm((eye(n) - B * B_inv) * (S * A - L * S), 'fro') < 1e-12
            K_final = B_inv * (A - inv(S) * L * S);
            eigs = eig(A - B * K_final);
            if all(abs(eigs) <= 1)
                flag = true;
                return;
            end
        end
        i = i + 1; 
        e(i) = min_norm;
    end
    
    if flag == false
        K_final = B_inv * (A - inv(S) * L * S);
        eigs = eig(A - B * K_final);
        if all(abs(eigs) <= 1)
            flag = true;
            return;
        end
        warning('Problem infeasible: no static feedback found by CVX.'); 
    end
end