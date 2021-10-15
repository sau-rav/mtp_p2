function L = projectNorm(L)
    val = norm(L, 'fro');
    delta = 1e-6;
    L = L / max(1 + delta, val + delta);
end