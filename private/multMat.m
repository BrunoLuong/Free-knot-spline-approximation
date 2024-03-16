function Balpha = multMat(B, alpha)
% %ake the matrix product B*alpha
% Balpha = B*alpha; without returning NaN

col = any(B,1);
Balpha = B(:,col)*alpha(col,:);
% Balpha = B*alpha;

end % multMat
