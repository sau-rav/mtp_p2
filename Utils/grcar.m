function A = grcar(n, k) 
    A = eye(n); 
    A = A - diag(ones(n - 1, 1), -1); 
    for i = 1 : k
        A = A + diag(ones(n - i, 1), i);
    end
end