function [S, L, K_final, flag] = testGradDesc(A, B, options)
    
    n = size(A, 1);
    B_inv = pinv(B);
    i = 1;
    K_final = [];
    flag = false;
    
    % Random Initialisation
    S = randn(n); S = S * S'; 
    L = rand(n); L = L / (5 * norm(L, 'fro'));
    
    min_norm = norm((eye(n) - B * B_inv) * (S * A - L * S), 'fro');
    fprintf('norm value initial : %2.10f\n', min_norm);
    e = [];
    e(i) = min_norm;
    
    disp('Solving BCD for checking feasibility');
    while i < options.itermax ...
        && (i <= 10 || e(i - 1) - e(i - 2) > 1e-6)
        
        fprintf('iteration : %d\n', i);
        % keep S as constant solve for L
        SSt = S * S';
        c = eye(n) - B * B_inv;
        
        LC1 = eigs(c, 1) * eigs(SSt, 1);
        step1 = 1 / LC1;
        gL = -norm(eye(n) - B * B_inv, 'fro')^2 * (S * A - L * S) * S';
        
        Ln = L - gL * step1; 
        Ln = projectNorm(Ln); 
        
        fprintf('norm value by subproblem 1 : %2.10f\n', norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro'));
        if norm((eye(n) - B * B_inv) * (S * A - Ln * S), 'fro') < e(i)
            fprintf('accepted value by subproblem 1\n');
            L = Ln;
        end
        
        % keep L as constant solve for S
        AAt = A * A';
        LtL = L' * L;
        
        LC2 = eigs(c, 1) * (eigs(AAt, 1) + eigs(LtL, 1) + 2 * eigs(L, 1) * eigs(A, 1));
        step2 = 1 / LC2;
        gS = norm(eye(n) - B * B_inv, 'fro')^2 * ((S * A - L * S) * A' - L' * (S * A - L * S));
        
        Sn = S - gS * step2;
        Sn = projectPSD(Sn);
        
        fprintf('norm value by subproblem 2 : %2.10f\n', norm((eye(n) - B * B_inv) * (Sn * A - L * Sn), 'fro'));
        if norm((eye(n) - B * B_inv) * (Sn * A - L * Sn), 'fro') < e(i)
            fprintf('accepted value by subproblem 2\n');
            S = Sn;
        end
        
        
        err = norm((eye(n) - B * B_inv) * (S * A - L * S), 'fro');
        min_norm = min(min_norm, err);
        fprintf('norm value at iteration end : %2.10f\n', err);
        
        if norm((eye(n) - B * B_inv) * (S * A - L * S), 'fro') < 1e-9
            K_final = B_inv * (A - inv(S) * L * S);
            evals = eig(A - B * K_final);
            if all(abs(evals) <= 1)
                flag = true;
                fprintf('norm value final : %2.10f\n', e(end));
                return;
            end
        end
        i = i + 1; 
        e(i) = err;
    end
    
    if flag == false
        K_final = B_inv * (A - inv(S) * L * S);
        evals = eig(A - B * K_final);
        if all(abs(evals) <= 1)
            flag = true;
            return;
        end
        warning('Problem infeasible: no static feedback found by gradient descend'); 
    end
    fprintf('norm value final : %2.10f\n', e(end));
end