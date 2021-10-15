% This file shows a numerical example of computing a solution to the 
% static-output feedback (SOF) problem with minimum norm, that is, to solve
%
% min_{K}  ||K||  such that A-BKC is stable for discrete time LTI system
% 
% We have used 5 initialization for SSDP method
% The initialization used are:
% 1. Random    - S = rand()
% 2. Standard  - S = eye() 
% 3. LMI       - S = sqrt(dlyap(A, eye(n)))
% 4. ABI       - S = solution of ABI
% 5. AIC       - S = solution of AIC
% 
% To display the graph : set figdisp to true
%% Initialisation of matrices
clear all;
clc;

addpath('./Utils/COMPlib_r1_1');
addpath('./Utils');
% testing feasibility
n = 4; m = 2; p = 2;

% A = rand(n);
% A = grcar(n, m);
% A = gallery('dramadah', n); 
% B = randn(n, m);
% C = randn(p, n);
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC1');

options.timemax = 20; 
options.itermax = 50;
% initialisation - 1 : S = eye(n), 2 : S = rand(), 3 : S = dlyap(), 4 : ABI, 5 : AIC 
options.init = 1;
options.display = 1;
fig = false;

%% feasibility test using SSDP with 5 different initialization
% result variables
feasRes = Inf * ones(1,5);
optRes = Inf * ones(1,5);
feasT = zeros(1,5);
optT = zeros(1,5);

for i = 1:5
     try
        options.init = i;
        fprintf('Checking feasibility of SOF for initialisation %d...\n', options.init);
        cput = cputime;
        [S, L, K, flag] = SOFFeas(A, B, C, options);
        feasT(i) = cputime - cput;
        feasRes(i) = norm(K, 'fro');
        fprintf('|| K || = %2.10f\n', feasRes(i));
        if options.display == 2
            fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
            fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
            
            disp(eig(A));
            disp(eig(A - B * K * C));
            fprintf('All eigen value of A - BKC inside unit circle : %d\n', all(abs(eig(A - B * K * C)) <= 1));
        end
        flag = flag & all(abs(eig(A - B * K * C)) <= 1);
        if fig == true
            figure 
            plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
            plot(eig(A - B * K * C) + 1e-16*sqrt(-1), 'ro'); hold on;
            angles = 0:2*pi/1000:2*pi;
            xa = cos(angles);
            ya = sin(angles);
            plot(xa, ya, 'k--'); 
            legend('eig(A)', 'eig(A - BKC)', 'Unit circle' ); 
        end
        
        if flag == false
            feasRes(i) = Inf;
            disp('Problem infeasible: No static feedback found');
        else  
            disp('Problem feasible : Static feedback found');
            cput = cputime;
            [S, L, K, ~, ~] = SOFMin(A, B, C, S, L, options);
            optT(i) = cputime - cput;
            optRes(i) = norm(K, 'fro');
            fprintf('|| K ||f = %2.10f\n', optRes(i));
            fprintf('All eigen value of A - BKC inside unit circle : %d\n', all(abs(eig(A - B * K * C)) <= 1));
            flag = flag & all(abs(eig(A - B * K * C)) <= 1);
            if flag == false
                optRes(i) = feasRes(i);
            end
            if fig == true
                figure 
                plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
                plot(eig(A - B * K * C) + 1e-16*sqrt(-1), 'ro'); hold on;
                angles = 0:2*pi/1000:2*pi;
                xa = cos(angles);
                ya = sin(angles);
                plot(xa, ya, 'k--'); 
                legend('eig(A)', 'eig(A - BKC)', 'Unit circle' );
            end
        end
     catch
        fprintf('Error for init %d\n', i);
        continue;
     end
end
fprintf('\n-------------------------------------------------------------\n');
fprintf('eye(n)  :   %2.3f\n', optRes(1));
fprintf('rand(n) :   %2.3f\n', optRes(2));
fprintf('dlyap() :   %2.3f\n', optRes(3));
fprintf('ABI     :   %2.3f\n', optRes(4));
fprintf('AIC     :   %2.3f\n', optRes(5));
fprintf('---------------------------------------------------------------\n');


% output for result section
% for i = 1:5
%     if(feasRes(i) == Inf) 
%         if 0.1 < feasRes(i) && 10.0 > feasRes(i)
%             fprintf('%1.2f ', feasRes(i));
%         else
%             fprintf('%1.2e ', feasRes(i));
%         end
%         if 0.1 < optRes(i) && 10.0 > optRes(i)
%             fprintf('%1.2f ', optRes(i));
%         else
%             fprintf('%1.2e ', optRes(i));
%         end
%     else
%         if 0.1 < feasRes(i) && 10.0 > feasRes(i)
%             fprintf('%1.2f(%2.2f) ', feasRes(i), feasT(i));
%         else
%             fprintf('%1.2e(%2.2f) ', feasRes(i), feasT(i));
%         end
%         if 0.1 < optRes(i) && 10.0 > optRes(i)
%             fprintf('%1.2f(%2.2f) ', optRes(i), optT(i));
%         else
%             fprintf('%1.2e(%2.2f) ', optRes(i), optT(i));
%         end
%     end
% end
% fprintf('\n');
% end
