function Balpha = multMat(B, alpha)
% Make the matrix "robust" product B*alpha
% Balpha = B*alpha;

if all(isfinite(alpha),'all')
    Balpha = B*alpha;
else
    % Balpha = B*alpha;
    % without returning NaN, since alpha ù*might contain Inf
    col = any(B,1);
    Balpha = B(:,col)*alpha(col,:);
end

if (issparse(B) || issparse(alpha)) % && ~issparse(Balpha)
    % This seems to help later, , not sure why (!), but otherwise
    % there might be a division by 0 in DerivBKnotDeriv()
    % testBSFK(3) sometime reveals the flaw
    Balpha = sparse(Balpha);
end

end % multMat

