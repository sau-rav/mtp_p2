function Qp = projectPSD(Q, epsilon, delta) 
    if isempty(Q)
        Qp = Q;
        return;
    end
    if nargin <= 1
        epsilon = 0;
    end
    if nargin <= 2
        delta = +Inf;
    end
    Q = (Q + Q') / 2; 
    if max(max(isnan(Q))) == 1 || max(max(isinf(Q))) == 1
        error('Input matrix has infinite or NaN entries');
    end
    [V, e] = eig(Q); 
    Qp = V * diag(min(delta, max(diag(e), epsilon))) * V'; 
end