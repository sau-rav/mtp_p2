% BCD + CVX method for solving the feasibility problem of 
% min_{K}  ||K||  such that A-BK is stable for discrete time LTI system
% 
% BCD method has been used to solve the using solutions to the convex
% subproblems created for :
% u = inf || (I - BB+) (AT - TL) ||
% where T = inv(S), T > 0 (positive definite) and || L || <= 1
%
% For solving the convex subproblem CVX has been used
%
% Input : A, B, options
% Output : S = inv(T) : positive definite matrix
%          L          : ||L|| <= 1
%          K          : (I - BB+)(A - inv(S)LS)
%          flag       : True if feasible False if infeasible
%          e          : error evolution over time
%          t          : timestamps
% ----------------------------------------------------------------------

function [S, L, K_final, flag, e, t] = SSFFeasCVX(A, B, options)
    n = size(A, 1);
    B_inv = pinv(B);
    projB  = (eye(n) - B * B_inv);
    K_final = Inf;
    i = 1;
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
    
    while i < options.itermax ...
        && (i <= 3 || e(i - 1) - e(i - 2) > 1e-6)
        
        e(i + 1) = e(i);
        % keep S as constant solve for L
        cvx_begin sdp quiet
            variable Ln(n,n) 
            minimize(norm(projB * (A * T - T * Ln), 'fro'))
            subject to
                norm(Ln, 'fro') <= 1;
        cvx_end
        if norm(projB * (A * T - T * Ln), 'fro') < e(i + 1)
            L = Ln;
            e(i + 1) = norm(projB * (A * T - T * L), 'fro');
        end
        
        % keep L as constand solve for S
        cvx_begin sdp quiet
            variable Tn(n,n)
            minimize(norm(projB * (A * Tn - Tn * L), 'fro'))
            subject to
                Tn == semidefinite(n); 
        cvx_end 
        if norm(projB * (A * Tn - Tn * L), 'fro') < e(i + 1)
            T = Tn;
            e(i + 1) = norm(projB * (A * T - T * L), 'fro');
        end
        
        if options.display == 2
            fprintf('norm value at iteration %d off %d : %2.10f\n', i, options.itermax, e(i + 1));
        end
        t(i + 1) = cputime - ts;
        S = pinv(T);
        if norm(projB * (A * T - T * L), 'fro') < 1e-12
            K_final = B_inv * (A - T * L * S);
            eig_vals = eig(A - B * K_final);
            if all(abs(eig_vals) <= 1)
                flag = true;
                return;
            end
        end
        i = i + 1; 
    end
    
    S = pinv(T);
    if flag == false
        K_final = B_inv * (A - T * L * S);
        eig_vals = eig(A - B * K_final);
        if all(abs(eig_vals) <= 1)
            flag = true;
            return;
        end
    end
end