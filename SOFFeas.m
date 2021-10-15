% SSDP method for checking for feasibility of 
% K = inf || (I - BB+) * (A - pinv(S) * L * S)|| + 
%                       || (A - pinv(S) * L * S) * (CC+ - I)|| 
% such that A-BKC is stable for discrete time LTI system
% 
% SSDP solves for Del(S) and Del(L) such that (S + Del(S), L + Del(L)) is
% a better solution than (S, L)
% using approximation final optimization problem becomes
% minimize norm( (I - BB+) * (A - (inv(S)LS + inv(S)LDel(S) + inv(S)Del(L)S 
%           - inv(S)Del(S)inv(S)Del(L)S) + (A - (inv(S)LS + inv(S)LDel(S) 
%           + inv(S)Del(L)S - inv(S)Del(S)inv(S)Del(L)S) * (CC+ - I)
%
% subject to
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

function [S, L, K_final, flag, err, t] = SOFFeas(A, B, C, options)
    n = size(A, 1);
    cput = cputime;
    B_inv = pinv(B); 
    projB  = (eye(n) - B * B_inv); 
    C_inv = pinv(C); 
    projC  = (eye(n) - C_inv * C);
    flag = false;
    
    %if A is stable 
    if all(abs(eigs(A)) <= 1)
        S = eye(n);
        L = A / (1.5 * norm(A, 'fro'));
        
        T = pinv(S);
        cvx_begin sdp quiet
            variable L(n,n)
            minimize(norm(projB * (A - T * L * S), 'fro') + ...
                     norm((A - T * L * S) * projC, 'fro'))
            subject to
                norm(L, 'fro') <= 1;
        cvx_end

        K_final = B_inv * (A - inv(S) * L * S) * C_inv;
        flag = true;
        return;
    end
    
	% initializations
    if options.init == 1
        S = eye(n);
        L = A / norm(A, 'fro');
    elseif options.init == 2
        S = randn(n); S = S * S'; 
        L = rand(n); L = L / (2 * norm(L));
    elseif options.init == 3
        P = dlyap(A, eye(n));
        S = sqrtm(P);
        singular_val = svd(A);
        L = A / singular_val(1);
    elseif options.init == 4
        optionsinit.itermax = 10; 
        optionsinit.init = 2;
        optionsinit.display = 1;
        [S, L, ~, ~, ~, ~] = SSFFeasCVX(A, B, optionsinit);
    elseif options.init == 5
        optionsinit.itermax = 10; 
        optionsinit.init = 2;
        optionsinit.display = 1;
        [S, L, ~, ~, ~, ~] = SSFFeasCVX(A', C', optionsinit);
    end
    
    % solve for L given S
    T = pinv(S);
    cvx_begin sdp quiet
        variable L(n,n)
        minimize(norm(projB * (A - T * L * S), 'fro') + ...
                 norm((A - T * L * S) * projC, 'fro'))
        subject to
            norm(L, 'fro') <= 1;
    cvx_end
    
    i = 1;
    err(i) = norm(projB * (A - T * L * S), 'fro') + norm((A - T * L * S) * projC, 'fro');
    t(i) = cputime - cput; 
    
    epsilon = 1;
    K_final = B_inv * (A - T * L * S) * C_inv; 
    while i <= options.itermax && err(i) > 1e-9 && ... 
             (i <= 2 || (err(i-2) - err(i-1)) > 1e-3 * err(i-1))
        
        err(i + 1) = +Inf; 
        while err(i+1) > err(i) && epsilon > 1e-9
            T = pinv(S);
            cvx_begin sdp quiet
                variable DL(n,n)
                variable DS(n,n)
                
                minimize(norm(projB * (A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)), 'fro') + ... 
                                norm((A - (T*L*S + T*L*DS + T*DL*S - T*DS*T*L*S)) * projC, 'fro')); 
                subject to
                    norm(L + DL, 'fro') <= 1; 
                    S + DS == semidefinite(n);
        
                    norm(DL, 'fro') <= epsilon * norm(L, 'fro'); 
                    norm(DS, 'fro') <= epsilon * norm(S, 'fro'); 
            cvx_end
                
            if norm(DS, 'fro') > 0
                Lnew = L + DL;
                Snew = S + DS;
                Tnew = pinv(Snew);
                err(i+1) = norm(projB * (A - Tnew * Lnew * Snew), 'fro') + ...
                    norm((A - Tnew * Lnew * Snew) * projC, 'fro');
                epsilon = epsilon/2; 
            else
                err(i+1) = Inf; 
                epsilon = epsilon/10; 
            end
            t(i+1) = cputime - cput; 
        end
        
        if epsilon < 1e-9
            K_final = B_inv * (A - pinv(S) * L * S) * C_inv; 
            if all(abs(eig(A - B * K_final * C)) <= 1)
                flag = true;
                return;
            else
				if options.display == 2
                	disp('The algorithm failed to stabilize (A,B,C).'); 
                end
				return;
            end
        end
        epsilon = min(1, epsilon*4);
        
        L = L + DL;
        S = S + DS; 
        K_final = B_inv * (A - pinv(S) * L * S) * C_inv; 
        if all(abs(eig(A - B * K_final * C)) <= 1)
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
