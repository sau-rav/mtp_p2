% BCD + Gradient Descent method for solving the feasibility problem of 
% min_{K}  ||K||  such that A-BK is stable for discrete time LTI system
% 
% BCD method has been used to solve the using solutions to the convex
% subproblems created for :
% u = inf || (I - BB+) (AT - TL) ||
% where T = inv(S), T > 0 (positive definite) and || L || <= 1
%
% For solving the convex subproblem Gradient Descent method has been used
%
% Input : A, B, options
% Output : S = inv(T) : positive definite matrix
%          L          : ||L|| <= 1
%          K_final    : (I - BB+)(A - inv(S)LS)
%          flag       : True if feasible False if infeasible
%          e          : error evolution over time
%          t          : timestamps
% ----------------------------------------------------------------------
function [S, L, K_final, flag, e, t] = SSFFeasGrad(A, B, options)
    
    n = size(A, 1);
    B_inv = pinv(B);
    projB  = (eye(n) - B * B_inv);
    i = 1;
    K_final = [];
    flag = false;
    
    %if A is stable 
    if all(abs(eigs(A)) <= 1)
        S = eye(n);
        L = A / (1.5 * norm(A, 'fro'));
        
        T = pinv(S);
        cvx_begin sdp quiet
            variable L(n,n)
            minimize (norm(projB * (A - T * L * S), 'fro'))
            subject to
                norm(L, 'fro') <= 1;
        cvx_end

        K_final = B_inv * (A - inv(S) * L * S);
        e = []; t = [];
        flag = true;
        return;
    end
    
    % Initialisation
    if options.init == 1
        S = rand(n); S = S * S'; 
        L = rand(n); L = L / (2 * norm(L));
    elseif options.init == 2
        S = eye(n);
        L = A / norm(A, 'fro');
    elseif options.init == 3
        P = dlyap(A, eye(n));
        S = sqrtm(P);
        singular_val = svd(A);
        L = A / singular_val(1);
    end
    
    T = pinv(S);
    e = []; t = [];
    e(i) = norm(projB * (A * T - T * L), 'fro');
    ts = cputime;
    t(i) = cputime - ts;
    if options.display == 2
        fprintf('norm value initial : %2.10f\n', e(i));
    end
    
    while i < options.itermax ...
        && (i <= 3 || e(i - 1) - e(i - 2) > 1e-6)
        
        % keep S as constant solve for L
        Ln = gradDescL(A, B, B_inv, T, L);
        if options.display == 2
            fprintf('norm value by subproblem 1 : %2.10f\n', norm(projB * (A * T - T * Ln), 'fro'));
        end
        if norm(projB * (A * T - T * Ln), 'fro') < e(i)
            if options.display == 2
                fprintf('accepted value by subproblem 1\n');
            end
            L = Ln;
        end
        
        % keep L as constant solve for S
        Tn = gradDescS(A, B, B_inv, T, L);
        if options.display == 2
            fprintf('norm value by subproblem 2 : %2.10f\n', norm(projB * (A * Tn - Tn * L), 'fro'));
        end
        if norm(projB * (A * Tn - Tn * L), 'fro') < e(i)
            if options.display == 2 
                fprintf('accepted value by subproblem 2\n');
            end
            T = Tn;
        end
        
        err = norm(projB * (A * T - T * L), 'fro');
        % min_norm = min(min_norm, err);
        if options.display == 2
            fprintf('norm value at iteration end : %2.10f\n', err);
        end
        
        t(i + 1) = cputime - ts;
        e(i + 1) = err;
        
        if norm(projB * (A * T - T * L), 'fro') < 1e-9
            K_final = B_inv * (A - T * L * pinv(T)); % inv(S) == T == pinv(S)
            eig_vals = eig(A - B * K_final);
            if all(abs(eig_vals) <= 1)
                flag = true;
                S = pinv(T);
                return;
            end
        end
        i = i + 1; 
    end
    
    if flag == false
        K_final = B_inv * (A - T * L * pinv(T));
        eig_vals = eig(A - B * K_final);
        if all(abs(eig_vals) <= 1)
            flag = true;
            S = pinv(T);
            return;
        end
    end
    S = pinv(T);
end