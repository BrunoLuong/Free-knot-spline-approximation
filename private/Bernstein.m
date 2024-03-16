function B = Bernstein(x, t, j, k, alpha, leftflag)
% B = Bernstein(x, t, j, k)
% B = Bernstein(x, t, j, k, alpha)
% Compute Bernstein polynomial basis using de De Boor's algorithm
%
% INPUTS:
%   x: vectors/array, point coordinates at which the function is to be
%      evaluated
%   t: vector, knots points, must be ascending sorted
%   j: vector, vector of spatial index, must be in [1:length(t)-k]
%      if it's empty all the basis functions are computed.
%   k-1: "order" of the spline (k is scalar)
%      k==1 -> piecewise constant
%      k==2 -> linear
%      k==4 -> cubic
%   alpha: vectors of size length(t)-k, optional coefficients of the basis
% OUTPUTS:
%   B: (m x n) where m is length(x), n is length(j)              
%       Each column of B is the basis function B_j,k
%   If t is not monotonically non-decreasing, B will be an empty array
%
%   Note: B_j,k has support on [t(j),t(j+k)[
%
%  Call with TRUE for 6th argument: B = Bernstein(..., TRUE) to compute
%  the "Left B-spline", which has support on ]t(j),t(j+k)].
%  The left B-spline is needed e.g. to evaluate recursion integral between
%  two B-splines.
%
% Author: Bruno Luong
%   31-May-2010: correct bug for LEFT Bspline is called with scalar
%   10-Jun-2010: Change in Bernstein.m to avoid NaN for end knots
%   13-Mar-2024: Cleaner code in recursion loop
%   14-Mar-2023: De Boor algorithm, recusrsion on unnormalized B-spline
%       basis function

%%
coefin = nargin>=5 && ~isempty(alpha);

if coefin
    if size(alpha,1)==1
        alpha = alpha.';
    end
end

cls = class(x);

%%
if isvector(x) %&& ~coefin
    szx = length(x);
else
    szx = size(x);
end
% sort x so that it fall into subinterval by contiguous segment
[x, isx] = sort(reshape(x, [], 1));
m = size(x,1);

%%
% Special case, we ignore the Dirac distribution
if k<=0    
    if coefin
        B = zeros([szx size(alpha,2)],cls);
    else
        B = zeros([szx numel(j)],cls);
    end
    return
end

%%
% Max possible value of indice
maxj = length(t)-k;
if isempty(j) || coefin
    % all the index j
    j = 1:maxj;
    js = j;
else
    js  = sort(j(:));
end

% left and right bracket
jmin = js(1);
jmax = js(end);
% Check
if jmin<1 || jmax>maxj
    error('BERNSTEIN: j must be within [%d,%d]', 1, maxj);
end

%% unnormalized B-Spline basis
% k=1: step functions (piecwise constant)
n = jmax+k-jmin;
M = zeros([m, n], cls);

if nargin>=6 && leftflag
    % Left B-spline
    tt = t(jmax+k:-1:jmin);
    if issorted(-tt)
        [~, col] = histc(-x,-tt); %#ok
        col = length(tt)-col; % Correct BUG, 31/05/2010
    else
        B = [];
        return
    end
else
    % Default, right B-spline
    tt = t(jmin:jmax+k);
    if issorted(tt)
        [~, col] = histc(x,tt); %#ok
    else
        B = [];
        return
    end
end

inside = col>=1 & col<=size(M,2);
i = col; % interval index, sams size as x
row = find(inside);
col = i(inside); % also i(row)

if k >= 2 && m % i must not be empty, required by GetBreaksArray

    % Construct the array of break indices of sorted x array
    breaks = GetBreaksArray(i, jmax+k);

    %  eqt (5), Carl de Boor "On Calculating with B-spline" 1972 paper
    dt = t(col+1)-t(col);
    idt = 1 ./ dt;
    idt(dt == 0) = 0;
    M(row+(col-1)*m) = idt;

    %% Loop on order
    for kk=2:k
        % sub-interval to be processed
        jv = jmin:jmax+k-kk;
        % support lengths
        dtv = t(jv+kk)-t(jv);
        % Cox–de Boor recursion formula
        % recursion loop on sub-intervals
        for c=1:size(jv,2)
            dt = dtv(c);
            if dt~=0
                jj = jv(c);
                % which points x fall inside the support
                r = breaks(jj):breaks(jj+kk)-1;
                w = (x(r)-t(jj)) / dt;
                % eqt (8), Carl de Boor "On Calculating with B-spline" 1972 paper
                M(r,c) = w.*M(r,c) + (1-w).*M(r,c+1);
            end
        end
    end

    % Normalized B-Spline basis, 
    % B here is noted by N in Carl de Boor "On Calculating with B-spline" 1972 paper
    % definition (4)
    B = zeros(size(M),cls);
    jv = jmin:jmax;
    dtv = t(jv+k)-t(jv);
    for c=1:size(jv,2)
        dt = dtv(c);
        if dt~=0
            jj = jv(c);
            % which points x fall inside the support
            r = breaks(jj):breaks(jj+k)-1;
            B(r,c) = M(r,c)*dt;
        end
    end
else
    B = zeros(size(M), cls);
    B(row+(col-1)*m) = 1;
end

if length(j) ~= size(B,2)
    % Map to original vector j
    %[tf loc] = ismemberc(j, jmin:jmax); %#ok
    loc = ismembc2(j, jmin:jmax);
    B = B(:,loc);
end

%% Restore the order of original x
if ~isequal(isx, (1:size(B,1)).')
    B = B(isx,:,:);
end

%%
if coefin
    % Compute function from the coefficients
    B = multMat(B, alpha); % Bug fix 10-Jun-2010
    B = reshape(B, [szx size(alpha,2)]);
else
    % Basis
    B = reshape(B, [szx numel(j)]);
end

end % Bernstein
