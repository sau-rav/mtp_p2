clear all;
clc;

addpath('./Utils/COMPlib_r1_1');
addpath('./Utils');
m = 6; n = 6;

% A = rand(m);
% A = grcar(m, n);
% A = gallery('dramadah', m);
% B = rand(m, n);
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC1');

options.timemax = 20; 
options.itermax = 10;
% initialisation - 1 : random, 2 : standard, 3 : LMI based 
options.init = 2;
% display - 1 : result, 2 : result and debug
options.display = 1;
figdisp = false;

cput = cputime;
[S, L, K, flag1, errF1, tF1] = SSFFeasSSDP(A, B, options);
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
    title('SSDP : Before optimisation');
end

if flag1 == false
    feas(4) = Inf;
    disp('Problem infeasible: No static feedback found by SSDP');
else  
    % solve for optimal K norm
    fprintf('Static feedback found by SSDP : Minimising the norm value...\n');
    options.algo = 2;
    cput = cputime;
    [S, L, K, errO4, tO4] = SSFMin(A, B, S, L, options);
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