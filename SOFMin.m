% SSDP method optimizing for 
% K = inf || inv(B) * (A - pinv(S) * L * S) * inv(S)||  
% such that A-BKC is stable for discrete time LTI system
% 
% SSDP solves for Del(S) and Del(L) such that (S + Del(S), L + Del(L)) is
% a better solution than (S, L)
% using approximation final optimization problem becomes
% minimize % SSDP method for checking for feasibility of 
% K = inf || (I - BB+) * (A - pinv(S) * L * S)|| + 
%                       || (A - pinv(S) * L * S) * (CC+ - I)|| 
% such that A-BKC is stable for discrete time LTI system
% 
% SSDP solves for Del(S) and Del(L) such that (S + Del(S), L + Del(L)) is
% a better solution than (S, L)
% using approximation final optimization problem becomes
% minimize norm( inv(B) * (A - (inv(S)LS + inv(S)LDel(S) + inv(S)Del(L)S 
%           - inv(S)Del(S)inv(S)Del(L)S)) * inv(C)
%
% subject to
% norm( (I - BB+) * (A - (inv(S)LS + inv(S)LDel(S) + inv(S)Del(L)S 
%           - inv(S)Del(S)inv(S)Del(L)S) + (A - (inv(S)LS + inv(S)LDel(S) 
%           + inv(S)Del(L)S - inv(S)Del(S)inv(S)Del(L)S) * (CC+ - I) < tol
% norm(L + Del(L)) <= 1; 
% S + Del(S) is positive definite; 
% norm(Del(S)) <= epsilon * norm(S);
% norm(Del(L)) <= epsilon * norm(L);
%
% Input : A, B, S, L, options
% Output : S          : positive definite matrix
%          L          : ||L|| <= 1
%          K_final    : inv(B)(A - inv(S)LS)inv(C)
%          flag       : true of feasible false if infeasible
%          err        : error evolution over time
%          t          : timestamps
% ----------------------------------------------------------------------

function [S, L, K_final, err, t] = SOFMin(A, B, C, S, L, options)
    
    n = size(A, 1);
    cput = cputime;
    B_inv = pinv(B); 
    projB  = (eye(n) - B * B_inv); 
    C_inv = pinv(C); 
    projC  = (eye(n) - C_inv * C);
    
    K_final = B_inv * (A - pinv(S) * L * S) * C_inv;
    tol1 = norm(projB * (A - pinv(S) * L * S), 'fro');
    tol2 = norm((A - pinv(S) * L * S) * projC, 'fro');
    i = 1;
    
    cput = cputime;
    err(i) = norm(K_final, 'fro');  
    t(i) = cputime-cput; 
    
    epsilon = 1;
 
    while i <= options.itermax && cputime - cput <= 2 * options.timemax && err(i) > 1e-9 ... 
            && (i <= 3 || (err(i-2)-err(i-1)) > 1e-3)
        err(i+1) = +Inf; 
        while  (isnan(err(i+1)) || err(i+1) > err(i)) && epsilon >= 1e-6
            T = pinv(S);
            DL = zeros(n, n);
            DS = zeros(n, n);
            
            try
                cvx_begin sdp quiet
                    variable DL(n,n)
                    variable DS(n,n)
                    minimize(norm(B_inv*(A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S))*C_inv, 'fro'))
                    
                    subject to
                    norm(projB * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol1;
                    norm((A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)) * projC, 'fro') <= tol2;
                    norm(DL, 'fro') <= epsilon * norm(L, 'fro');
                    norm(DS, 'fro') <= epsilon * norm(S, 'fro');
                    norm(L + DL, 'fro') <= 1; 
                    S + DS == semidefinite(n);
                cvx_end
            catch
                if options.display == 2
                    fprintf('Norm value cant be minimised further\n');
                end
                cvx_clear;
                return;
            end
            
            if (norm(projB*(A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') <= tol1 && ...
                    norm((A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S))*projC, 'fro') <= tol2) && ...
                    all(abs(eigs(A - B * (B_inv * (A - pinv(S + DS) * (L + DL) * (S + DS)) * C_inv) * C)) <= 1)
                Lnew = (L + DL);
                Snew = (S + DS);
                err(i+1) = norm(B_inv * (A - pinv(Snew) * Lnew * Snew) * C_inv, 'fro');
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
            K = B_inv * (A - pinv(S) * L * S) * C_inv;
            return;
        end
		if options.display == 2
        	fprintf('error after iteration %d : %2.10f\n', i + 1, err(i + 1));
		end
        K_final = B_inv * (A - pinv(S) * L * S) * C_inv; 
        t(i + 1) = cputime - cput;
        i = i+1; 
    end
end
