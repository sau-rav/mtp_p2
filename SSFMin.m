% SSDP method for solving the optimization of 
% K = inf || B+ * (A - pinv(S) * L * S) * C+|| such that A-BK is stable 
% for discrete time LTI system
% 
% SSDP solves for Del(S) and Del(L) such that (S + Del(S), L + Del(L)) is
% a better solution than (S, L)
% using approximation final optimization problem becomes
% minimize norm( B_inv * (A - (inv(S)LS + inv(S)LDel(S) + inv(S)Del(L)S 
%                                          - inv(S)Del(S)inv(S)Del(L)S))
                        
% subject to
% norm((eye(n) - B*B_inv) * (A - (inv(S)LS + inv(S)LDel(S) + inv(S)Del(L)S 
%                             - inv(S)Del(S)inv(S)Del(L)S))) <= tol; 
% norm(L + Del(L)) <= 1; 
% S + Del(S) is positive definite; 
% norm(Del(S)) <= epsilon * norm(S);
% norm(Del(L)) <= epsilon * norm(L);
%
% Input : A, B, S, L, options
% Output : S          : positive definite matrix
%          L          : ||L|| <= 1
%          K_final    : (B+)(A - inv(S)LS)
%          err        : error evolution over time
%          t          : timestamps
% ----------------------------------------------------------------------

function [S, L, K_final, err, t] = SSFMin(A, B, S, L, options)
    
    n = size(A, 1);
    B_inv = pinv(B);
    i = 1;
    K_final = B_inv * (A - pinv(S) * L * S);
    T = pinv(S);
    tol = norm((eye(n) - B * B_inv) * (A * T - T * L), 'fro');
    mult = 10;

    cput = cputime;
    err(i) = norm(K_final, 'fro');  
    t(i) = cputime-cput; 
    
    epsilon = 1;
 
    while i <= options.itermax && cputime - cput <= 2 * options.timemax && err(i) > 1e-9 ... 
            && (i <= 3 || (err(i-2)-err(i-1)) > 1e-3)
        err(i+1) = +Inf; 
        while  (isnan(err(i+1)) || err(i+1) > err(i)) && epsilon >= 1e-6
            if rem(i, 2) == 0
                T = pinv(S);
                DL = zeros(n, n);
                DS = zeros(n, n);
                
                try
                    cvx_begin sdp quiet
                        variable DS(n,n)
                        minimize(norm(B_inv*(A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro'))
                        
                        subject to
                        norm((eye(n) - B*B_inv) * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol; 
                        S + DS == semidefinite(n); 
                        norm(DS, 'fro') <= epsilon * norm(S, 'fro');
                    cvx_end
                                        
                    cvx_begin sdp quiet
                        variable DL(n,n)
                        minimize(norm(B_inv*(A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro'))
                        
                        subject to
                        norm((eye(n) - B*B_inv) * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol; 
                        norm(L + DL, 'fro') <= 1; 
                        norm(DL, 'fro') <= epsilon * norm(L, 'fro');
                    cvx_end 
                catch
                    if options.display == 2
                        fprintf('Norm value cant be minimised further\n');
                    end
                    cvx_clear;
                    return;
                end
            else
                T = pinv(S);
                DL = zeros(n, n);
                DS = zeros(n, n);
                try
                    cvx_begin sdp quiet
                        variable DL(n,n)
                        variable DS(n,n)
                        
                        minimize(norm(B_inv*(A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro'))
                        
                        subject to
                            norm((eye(n) - B*B_inv) * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol; 
                            norm(L + DL, 'fro') <= 1; 
                            S + DS == semidefinite(n); 
                            
                            norm(DS, 'fro') <= epsilon * norm(S, 'fro');
                            norm(DL, 'fro') <= epsilon * norm(L, 'fro');
                    cvx_end
                catch
                    if options.display == 2
                        fprintf('Norm value cant be minimised further\n');
                    end
                    cvx_clear;
                    return;
                end
            end
            
            if norm((eye(n) - B*B_inv) * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol*mult && ...
                    all(abs(eigs(A - B * (B_inv * (A - pinv(S + DS) * (L + DL) * (S + DS))))) <= 1)
                Lnew = (L + DL);
                Snew = (S + DS);
                err(i+1) = norm(B_inv * (A - pinv(Snew) * Lnew * Snew), 'fro');
                t(i+1) = cputime - cput;
                epsilon = epsilon/2; 
            else
                err(i+1) = +Inf;
                t(i+1) = cputime-cput;
                epsilon = epsilon/10; 
            end
        end
        
        epsilon = min(1, epsilon * 4); 
        if err(i+1) < err(i) 
            L = L + DL;
            S = S + DS;
        else
            if options.display == 2
                disp('The algorithm failed to return further descent direction'); 
            end
            K = B_inv * (A - pinv(S) * L * S);
            return;
        end
        if options.display == 2
            fprintf('error after iteration %d : %2.10f\n', i + 1, err(i + 1));
        end
        K_final = B_inv * (A - pinv(S) * L * S); 
        t(i + 1) = cputime - cput;
        i = i+1; 
    end
end
