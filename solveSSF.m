% This file shows a numerical example of computing a solution to the 
% static-state feedback (SSF) problem with minimum norm, that is, to solve
%
% min_{K}  ||K||  such that A-BK is stable for discrete time LTI system
% 
% We have used 4 methods with 3 initialization each
% The methods used are
% 1. BCD + CVX
% 2. BCD + Gradient Descent
% 3. BCD + Fast Gradient Descent
% 4. SSDP Method
%
% The initialization used are:
% 1. Random    - S = rand(), L = rand()
% 2. Standard  - S = eye() , L = A / (2 * norm(L))
% 3. LMI       - S = sqrt(dlyap(A, eye(n))), L = A / svd(A, 1)
% 
% To display the graph : set figdisp to true
%% Initialisation of matrices
clear all;
clc;

addpath('./Utils/COMPlib_r1_1');
addpath('./Utils');
m = 10; n = 6;

% A = rand(m);
A = grcar(m, n);
% A = gallery('dramadah', m);
B = rand(m, n);
% [A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC1');

options.timemax = 20; 
options.itermax = 10;
% initialisation - 1 : random, 2 : standard, 3 : LMI based 
options.init = 2;
% display - 1 : result, 2 : result and debug
options.display = 1;
figdisp = false; % change to true for figure display at end

feasRes = Inf * ones(1, 12);
optRes = Inf * ones(1, 12);
optT = zeros(1, 12);
feasT = zeros(1, 12);
% iterating over the initializations
for i = 1:3
    feas = Inf * ones(1, 4);
    opt = Inf * ones(1, 4);
    feasTime = zeros(1, 4);
    optTime = zeros(1, 4);
    options.init = i;
    
    % feasibility test using BCD + CVX
    fprintf('Feasibility problem using BCD + CVX ...\n');
    try
        cput = cputime;
        [S, L, K, flag1, errF1, tF1] = SSFFeasCVX(A, B, options);
        feasTime(1) = cputime - cput;
        feas(1) = norm(K, 'fro');
        fprintf('|| K || = %2.10f\n', feas(1));
        
        if options.display == 2
            fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
            fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
            fprintf('eig(A)'); disp(eig(A));
            fprintf('eig(A - BK)'); disp(eig(A - B * K));
            fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
        end
        flag1 = flag1 & all(abs(eig(A - B * K)) <= 1);
        if figdisp 
            figure 
            plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
            plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
            angles = 0:2*pi/1000:2*pi;
            xa = cos(angles);
            ya = sin(angles);
            plot(xa, ya, 'k--'); 
            legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
            title('CVX : Before optimisation');
        end
        
        if flag1 == false
            feas(1) = Inf;
            disp('Problem infeasible: No static feedback found by CVX');
        else  
            % solve for optimal K norm
            fprintf('Static feedback found by CVX : Minimising the norm value...\n');
            options.algo = 1;
            cput = cputime;
            [S, L, K, errO1, tO1] = SSFMin(A, B, S, L, options);
            optTime(1) = cputime - cput;
            if all(abs(eig(A - B * K)) <= 1)
                opt(1) = norm(K, 'fro');
            end
            fprintf('|| K || = %2.10f\n', opt(1));
            if options.display == 2
                fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
                fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
                fprintf('eig(A - BK)'); disp(eig(A - B * K));
                fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
            end
            if figdisp
                figure 
                plot(eig(A) + 1e-16*sqrt(-1), 'bo'); hold on; 
                plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
                angles = 0:2*pi/1000:2*pi;
                xa = cos(angles);
                ya = sin(angles);
                plot(xa, ya, 'k--'); 
                legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
                title('CVX : Post optimisation');
            end
        end
    catch
        fprintf('Error encountered for CVX : initialisation %d\n', i);
    end
    
    % feasibility test using BCD + Gradient Descent
    fprintf('\nFeasibility problem using BCD + Grad. Descent ...\n');
    try
        cput = cputime;
        [S, L, K, flag2, errF2, tF2] = SSFFeasGrad(A, B, options);
        feasTime(2) = cputime - cput;
        feas(2) = norm(K, 'fro');
        fprintf('|| K || = %2.10f\n', feas(2));
        
        if options.display == 2
            fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
            fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
            fprintf('eig(A)'); disp(eig(A));
            fprintf('eig(A - BK)'); disp(eig(A - B * K));
            fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
        end
        flag2 = flag2 & all(abs(eig(A - B * K)) <= 1);
        if figdisp
            figure 
            plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
            plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
            angles = 0:2*pi/1000:2*pi;
            xa = cos(angles);
            ya = sin(angles);
            plot(xa, ya, 'k--'); 
            legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
            title('Grad. Desc : Before optimisation');
        end
        
        if flag2 == false
            feas(2) = Inf;
            disp('Problem infeasible: No static feedback found by Gradient Descent');
        else  
            % solve for optimal K norm
            fprintf('Static feedback found by Gradient Descent : Minimising the norm value...\n');
            options.algo = 1;
            cput = cputime;
            [S, L, K, errO2, tO2] = SSFMin(A, B, S, L, options);
            optTime(2) = cputime - cput;
            if all(abs(eig(A - B * K)) <= 1)
                opt(2) = norm(K, 'fro');
            end
            fprintf('|| K || = %2.10f\n', opt(2));
            if options.display == 2
                fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
                fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
                fprintf('eig(A - BK)'); disp(eig(A - B * K));
                fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
            end
            if figdisp
                figure 
                plot(eig(A) + 1e-16*sqrt(-1), 'bo'); hold on; 
                plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
                angles = 0:2*pi/1000:2*pi;
                xa = cos(angles);
                ya = sin(angles);
                plot(xa, ya, 'k--'); 
                legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
                title('Grad. Desc : Post optimisation');
            end
        end
    catch
        fprintf('Error encountered for grad : initialisation %d\n', i);
    end
    
    % feasibility test using BCD + Fast Gradient Descent
    fprintf('\nFeasibility problem using BCD + Fast Grad. Descent ...\n');
    try
        cput = cputime;
        [S, L, K, flag3, errF3, tF3] = SSFFeasFast(A, B, options);
        feasTime(3) = cputime - cput;
        feas(3) = norm(K, 'fro');
        fprintf('|| K || = %2.10f\n', feas(3));
        
        if options.display == 2
            fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
            fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
            fprintf('eig(A)'); disp(eig(A));
            fprintf('eig(A - BK)'); disp(eig(A - B * K));
            fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
        end
        flag3 = flag3 & all(abs(eig(A - B * K)) <= 1);
        if figdisp
            figure 
            plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
            plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
            angles = 0:2*pi/1000:2*pi;
            xa = cos(angles);
            ya = sin(angles);
            plot(xa, ya, 'k--'); 
            legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
            title('Fast Grad. Desc : Before optimisation');
        end
        if flag3 == false
            feas(3) = Inf;
            disp('Problem infeasible: No static feedback found by Fast Grad. Descent');
        else  
            % solve for optimal K norm
            fprintf('Static feedback found by Fast Grad. Descent : Minimising the norm value...\n');
            options.algo = 1;
            cput = cputime;
            [S, L, K, errO3, tO3] = SSFMin(A, B, S, L, options);
            optTime(3) = cputime - cput;
            if all(abs(eig(A - B * K)) <= 1)
                opt(3) = norm(K, 'fro');
            end
            fprintf('|| K || = %2.10f\n', opt(3));
            if options.display == 2
                fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
                fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
                fprintf('eig(A - BK)'); disp(eig(A - B * K));
                fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
            end
            if figdisp
                figure 
                plot(eig(A) + 1e-16*sqrt(-1), 'bo'); hold on; 
                plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
                angles = 0:2*pi/1000:2*pi;
                xa = cos(angles);
                ya = sin(angles);
                plot(xa, ya, 'k--'); 
                legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
                title('Grad. Desc : Post optimisation');
            end
        end
    catch
        fprintf('Error encountered for fast grad : initialisation %d\n', i);
    end
    
    % feasibility test using SSDP
    fprintf('Feasibility problem using SSDP ...\n');
    try
        cput = cputime;
        options.init = i;
        [S, L, K, flag4, errF4, tO4] = SSFFeasSSDP(A, B, options);
        feasTime(4) = cputime - cput;
        feas(4) = norm(K, 'fro');
        fprintf('|| K || = %2.10f\n', feas(4));
        
        if options.display == 2
            fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
            fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
            fprintf('eig(A)'); disp(eig(A));
            fprintf('eig(A - BK)'); disp(eig(A - B * K));
            fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
        end
        flag4 = flag4 & all(abs(eig(A - B * K)) <= 1);
        if figdisp 
            figure 
            plot(eig(A) + 1e-16*sqrt(-1),'bo'); hold on; 
            plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
            angles = 0:2*pi/1000:2*pi;
            xa = cos(angles);
            ya = sin(angles);
            plot(xa, ya, 'k--'); 
            legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
            title('SSDP : Before optimisation');
        end
        
        if flag4 == false
            feas(4) = Inf;
            disp('Problem infeasible: No static feedback found by SSDP');
        else  
            % solve for optimal K norm
            fprintf('Static feedback found by SSDP : Minimising the norm value...\n');
            options.algo = 2;
            cput = cputime;
            [S, L, K, errF4, tO4] = SSFMin(A, B, S, L, options);
            optTime(4) = cputime - cput;
            if all(abs(eig(A - B * K)) <= 1)
                opt(4) = norm(K, 'fro');
            end
            fprintf('|| K || = %2.10f\n', opt(4));
            if options.display == 2
                fprintf('Positive definite S : %d\n', all(abs(eig(S)) > 0));
                fprintf('Norm of L <= 1 : %d\n', norm(L, 'fro') <= 1);
                fprintf('eig(A - BK)'); disp(eig(A - B * K));
                fprintf('All eigen value of A - BK inside unit circle : %d\n', all(abs(eig(A - B * K)) <= 1));
            end
            if figdisp
                figure 
                plot(eig(A) + 1e-16*sqrt(-1), 'bo'); hold on; 
                plot(eig(A - B * K) + 1e-16*sqrt(-1), 'ro'); hold on;
                angles = 0:2*pi/1000:2*pi;
                xa = cos(angles);
                ya = sin(angles);
                plot(xa, ya, 'k--'); 
                legend('eig(A)', 'eig(A - BK)', 'Unit circle' ); 
                title('SSDP : Post optimisation');
            end
        end
    catch 
        fprintf('Error encountered for SSDP : initialisation %d\n', i);
    end
    
    feasRes(i) = feas(1); feasT(i) = feasTime(1);
    feasRes(3 + i) = feas(2); feasT(3 + i) = feasTime(2);
    feasRes(6 + i) = feas(3); feasT(6 + i) = feasTime(3);
    feasRes(9 + i) = feas(4); feasT(9 + i) = feasTime(4);
    optRes(i) = opt(1); optT(i) = optTime(1);
    optRes(3 + i) = opt(2); optT(3 + i) = optTime(2);
    optRes(6 + i) = opt(3); optT(6 + i) = optTime(3);
    optRes(9 + i) = opt(4); optT(9 + i) = optTime(4);
    fprintf('\n-------------------------------------------------------------\n');
    fprintf('CVX       :   %2.3f\n', opt(1));
    fprintf('Grad      :   %2.3f\n', opt(2));
    fprintf('Fast Grad :   %2.3f\n', opt(3));
    fprintf('SSDP      :   %2.3f\n', opt(4));
    fprintf('---------------------------------------------------------------\n');

end


% used for printing results
% fprintf('CVX\n');
% for i = 1:3, fprintf('%2.5f ', feasRes(i)); end
% fprintf('\n');
% for i = 1:3, fprintf('%2.5f ', optRes(i)); end
% fprintf('\nGrad\n');
% for i = 4:6, fprintf('%2.5f ', feasRes(i)); end
% fprintf('\n');
% for i = 4:6, fprintf('%2.5f ', optRes(i)); end
% fprintf('\nFast Grad\n');
% for i = 7:9, fprintf('%2.5f ', feasRes(i)); end
% fprintf('\n');
% for i = 7:9, fprintf('%2.5f ', optRes(i)); end
% fprintf('\nSSDP\n');
% for i = 10:12, fprintf('%2.5f ', feasRes(i)); end
% fprintf('\n');
% for i = 10:12, fprintf('%2.5f ', optRes(i)); end
% fprintf('\n');


for i = 1:12
    if(feasRes(i) == Inf) 
        if 0.1 < feasRes(i) && 10.0 > feasRes(i)
            fprintf('%1.2f ', feasRes(i));
        else
            fprintf('%1.2e ', feasRes(i));
        end
        if 0.1 < optRes(i) && 10.0 > optRes(i)
            fprintf('%1.2f ', optRes(i));
        else
            fprintf('%1.2e ', optRes(i));
        end
    else
        if 0.1 < feasRes(i) && 10.0 > feasRes(i)
            fprintf('%1.2f(%2.2f) ', feasRes(i), feasT(i));
        else
            fprintf('%1.2e(%2.2f) ', feasRes(i), feasT(i));
        end
        if 0.1 < optRes(i) && 10.0 > optRes(i)
            fprintf('%1.2f(%2.2f) ', optRes(i), optT(i));
        else
            fprintf('%1.2e(%2.2f) ', optRes(i), optT(i));
        end
    end
end
fprintf('\n');
Co = ctrb(A,B);
unco = length(A) - rank(Co)