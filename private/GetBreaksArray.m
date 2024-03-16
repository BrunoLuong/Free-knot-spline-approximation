function breaks = GetBreaksArray(i, nknots)
% Construct the array of break indices of sorted x array
% INPUTS:
%   i: column vector of interval indices, must not be empty and sorted
%   nknots: scalar number of knots
% OUTUT:
%   breaks: column vector of length nknots
%       breaks(j) contains first index posistion in i such that
%       i(breaks(j)) >= j, for j=1,2,...,nknots

jumpidx = find([true; diff(i)~=0; true]);
i(end+1) = nknots+1; % wanted this value for last element of ij
                     % +1 so we be sure it will not repeat i(end)
ij = i(jumpidx);
breaks = interp1(ij, jumpidx, 1:nknots, 'next', 'extrap');

end % GetBreaksArray
