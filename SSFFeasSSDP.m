% SSDP method for solving the feasibility problem of 
% min_{K}  ||K||  such that A-BK is stable for discrete time LTI system
% 
% SSDP solves for Del(S) and Del(L) such that (S + Del(S), L + Del(L)) is
% a better solution than (S, L)
% u = inf || (I - BB+) (AT - TL) ||
% where T = inv(S), T > 0 (positive definite) and || L || <= 1
%
% Input : A, B, options
% Output : S = inv(T) : positive definite matrix
%          L          : ||L|| <= 1
%          K_final    : (I - BB+)(A - inv(S)LS)
%          flag       : True if feasible False if infeasible
%          err        : error evolution over time
%          t          : timestamps
% ----------------------------------------------------------------------

function [S, L, K_final, flag, err, t] = SSFFeasSSDP(A, B, options)
    n = size(A, 1);
    B_inv = pinv(B);
    projB  = (eye(n) - B * B_inv);
    K_final = Inf;
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

        K_final = B_inv * (A - pinv(S) * L * S);
        err = []; t = [];
        flag = true;
        return;
    end
    
    % initialisaiton for first approximate solution
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
    
    % solve for L given S
    % SSDP requires and initial solution so we solve for L keeping S 
    % constant using CVX
    T = pinv(S);
    cvx_begin sdp quiet
        variable L(n,n)
        minimize(norm(projB * (A - T * L * S), 'fro'))
        subject to
            norm(L, 'fro') <= 1;
    cvx_end
    
    i = 1;
    err = []; t = [];
    err(i) = norm(projB * (A - T * L * S), 'fro');
    cput = cputime;
    t(i) = cputime - cput; 
    
    epsilon = 1;
    K_final = B_inv * (A - T * L * S); 
    while i <= options.itermax && err(i) > 1e-9 && ... 
             (i <= 2 || (err(i-2) - err(i-1)) > 1e-3 * err(i-1))
        
        err(i + 1) = +Inf; 
        while err(i+1) > err(i) && epsilon > 1e-9
            T = pinv(S);
            cvx_begin sdp quiet
                variable DL(n,n)
                variable DS(n,n)
                
                minimize(norm(projB * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro')) 
                subject to
                    norm(L + DL, 'fro') <= 1; 
                    S + DS == semidefinite(n);
        
                    norm(DL, 'fro') <= epsilon * norm(L, 'fro'); 
                    norm(DS, 'fro') <= epsilon * norm(S, 'fro'); 
            cvx_end
                
            if norm(projB * (A - pinv(S+DS)*(L+DL)*(S+DS)), 'fro') < ...
                    norm(projB * (A - T * L * S), 'fro')
                Lnew = L + DL;
                Snew = S + DS;
                Tnew = pinv(Snew);
                err(i+1) = norm(projB * (A - Tnew * Lnew * Snew), 'fro');
                epsilon = epsilon/2; 
            else
                err(i+1) = Inf; 
                epsilon = epsilon/10; 
            end
            t(i+1) = cputime - cput; 
        end
        
        if epsilon < 1e-9
            K_final = B_inv * (A - pinv(S) * L * S); 
            if all(abs(eig(A - B * K_final)) <= 1)
                flag = true;
                return;
            else
                disp('The algorithm failed to stabilize (A,B).'); 
                return;
            end
        end
        epsilon = min(1, epsilon*4);
        
        L = L + DL;
        S = S + DS;  
        K_final = B_inv * (A - pinv(S) * L * S); 
        if all(abs(eig(A - B * K_final)) <= 1)
            flag = true;
            return;
        end
        if options.display == 2
        	fprintf('norm value at iteration %d off %d : %2.10f\n', i, options.itermax, err(i + 1));
		end
        i = i + 1;
    end
    flag = false;
end
