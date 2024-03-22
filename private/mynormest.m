function nA = mynormest(A)
% nA = mynormest(A)
% Estimate 2-norm of matrix A

if issparse(A)
    nA = normest(A);
else
    nA = norm(A);
end
end