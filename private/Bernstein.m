function B = Bernstein(x, t, j, k, alpha, leftflag)
% B = Bernstein(x, t, j, k)
% B = Bernstein(x, t, j, k, alpha)
% Compute Bernstein polynomial basis using de De Boor's algorithm
%
% INPUTS:
%   x: vectors/array, point coordinates at which the function is to be
%      evaluated
%   t: vector, knots positions, must be ascending sorted
%       usually t has the first and last elements repeat k times.
%   j: vector, vector of spatial index, must be in [1:length(t)-k]
%      if it's empty all the basis functions are computed.
%   k: scalar >= 1, "order" of the spline
%      k-1 is the degree of the polynomials on subintervals
%      k==1 -> piecewise constant
%      k==2 -> linear
%      k==4 -> cubic
%   alpha: vectors of length n := length(t)-k, optional coefficients of the
%       basis or matrix of the size (n x ndim). It can be seen as n spline
%       control points in a space of R^ndim
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
%   14-Mar-2024: De Boor algorithm, recusrsion on unnormalized B-spline
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
    js  = min(j):max(j);
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
    tt = -t(jmax+k:-1:jmin);
    if issorted(tt)
        col = length(tt)-myhistc(-x, tt);
    else
        B = [];
        return
    end
else
    % Default, right B-spline
    tt = t(jmin:jmax+k);
    if issorted(tt)
        col = myhistc(x, tt);
    else
        B = [];
        return
    end
end

% Construct the array of break indices of sorted x array
breaks = GetBreaksArray(col, length(tt));
GetSegment = @(j,l)GetSegmentHelper(j,l,breaks,jmin);

if k >= 2

    inside = col>=1 & col<=n;
    row = find(inside);
    col = col(inside); % also col(row)

    %  eqt (5), Carl de Boor "On Calculating with B-spline" 1972 paper
    c1 = jmin+col;
    c0 = c1-1;
    dt = t(c1)-t(c0);
    idt = 1 ./ dt;
    idt(dt == 0) = 0;
    M(row+(col-1)*m) = idt;

    %% Loop on order
    for kk=2:k
        % sub-interval to be processed
        jv = jmin:jmax+k-kk; % relative to full knot t
        % support lengths
        dtv = t(jv+kk)-t(jv);
        % Cox–de Boor recursion formula
        % recursion loop on sub-intervals
        for c=1:size(jv,2)
            dt = dtv(c);
            if dt %~=0
                jj = jv(c);
                % which points x fall inside the support
                r = GetSegment(jj,kk);
                w = (x(r)-t(jj)) / dt;
                % eqt (8), Carl de Boor "On Calculating with B-spline" 1972 paper
                M(r,c) = w.*M(r,c) + (1-w).*M(r,c+1);
            end
        end
    end

    % Normalized B-Spline basis, 
    % B here is noted by N in Carl de Boor "On Calculating with B-spline" 1972 paper
    % definition (4)
    n = length(j);
    B = zeros([m n],cls);
    dtv = t(j+k)-t(j);
    for c=1:n
        dt = dtv(c);
        if dt %~=0
            jj = j(c);
            cM = jj-jmin+1;
            % which points x fall inside the support
            r = GetSegment(jj,k);
            B(r,c) = M(r,cM)*dt;
        end
    end
else % k == 1
    n = length(j);
    B = zeros([m n],cls);
    for c=1:n
        jj = j(c);
        % which points x fall inside the support
        r = GetSegment(jj,k);
        B(r,c) = 1;
    end
end

%% Restore the order of original x
if ~isequal(isx, (1:size(B,1)).')
    B = B(isx,:,:);
end

%%
if coefin
    % Compute function from the coefficients
    B = multMat(B, alpha); % Bug fix 10-Jun-2010 no longer applied?
    B = reshape(B, [szx size(alpha,2)]);
else
    % Basis
    B = reshape(B, [szx numel(j)]);
end

end % Bernstein
