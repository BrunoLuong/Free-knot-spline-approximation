function [B, dB, dalpha] = BernKnotDeriv(x, knots, j, k, dknots, alpha, lambda, ...
                                       leftflag) %#codegen
% [B, dB] = BernKnotDeriv(x, knots, j, k, dknots)
% [s, ds] = BernKnotDeriv(x, knots, j, k, dknots, alpha)
% [s, ds, dalpha] = BernKnotDeriv(x, knots, j, k, dknots, alpha, lambda)
%
% Compute Bernstein polynomial basis by de De Casteljau's algorithm
% and its derivative with respect to the knots
%
% INPUTS:
%   x: vectors/array, point coordinates at which the function is to be
%      evaluated
%   knots: vector, knots points, must be ascending sorted
%   j: vector, vector of spatial index, must be in [1:length(knots)-k]
%         if it's empty all the basis functions are computed.
%   k-1: "order" of the spline (k is scalar)
%      k==1 -> piecewise constant
%      k==2 -> linear
%      k==4 -> cubic
%   dknots: increment of knots, column vector or 
%       array, each column is the direction where the derivative is computed.
%       In order to compute the Jacobian user must provide a basis vectors
%       for DKNOTS
%   alpha: vectors of size n:=length(knots)-k, optional coefficients of the
%          spline basis
%   lambda: dual variable, same number of elements as x
%
% OUTPUTS:
%   B: (m x n) where m is length(x), n is length(j)              
%       Each column of B is the basis function B_j,k
%   dB: the derivative of B with respect to the direction given by
%       dknots
%   s: spline function, same dimension as x
%   ds: the derivative of s with respcte to knots
%   dalpha: left product the derivative dB (of B) with the dual lambda,
%           dimension (n x p)
%
%   If knots is not monotonically non-decreasing, output will be an empty arrays
%
%   Note: B_j,k has support on [knots(j),knots(j+k)[
%
%  Call with TRUE for 8th argument: B = BernKnotDeriv(..., TRUE) to compute
%  the "Left B-spline", which has support on ]knots(j),knots(j+k)].
%  The left B-spline is needed e.g. to evaluate recursion integral between
%  two B-splines.      
%
% Author: Bruno Luong
%   31-May-2010: correct bug for LEFT Bspline is called with scalar
%   16-Mar-2023: De Boor algorithm, recusrsion on unnormalized B-spline
%       basis function

% Check if the coefficients are provided by user
coefin = nargin>=6 && ~isempty(alpha);

if coefin
    if size(alpha,1)==1
        alpha = alpha.';
    end
end
    
% Class of x
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
% Max possible value of indice
maxj = length(knots)-k;
if isempty(j) || coefin
    % all the index j
    j = 1:maxj;
    js = j;
else
    js  = sort(j(:));
end

% Reshape in column vector
if isequal(size(dknots),[1 length(knots)])
    dknots = dknots(:);
end
if isempty(dknots)
    dknots = zeros(length(knots),0,cls);
end
p = size(dknots,2);
dknots = reshape(dknots,[1 size(dknots)]);

% left and right bracket
jmin = js(1);
jmax = js(end);
% Check
if jmin<1 || jmax>maxj
    error('BERNSTEIN: j must be within [%d,%d]', 1, maxj);
end

%%
% Spcial case, we ignore the Dirac distribution
if k<=0    
    if coefin     
        B = zeros([szx size(alpha,2)],cls);
        dB = zeros([szx size(alpha,2) p],cls);
    else
        B = zeros([szx numel(j)],cls);
        dB = zeros([szx numel(j) p],cls);
    end
    return
end

%% unnormalized B-Spline basis
% k=1: step functions (piecwise constant)
n = jmax+k-jmin;
M = zeros([m, n], cls);

if nargin>=8 && leftflag
    % Left B-spline
    tt = knots(jmax+k:-1:jmin);
    if issorted(-tt)
        [~, col] = histc(-x,-tt); %#ok
        col = length(tt)-col; % Correct BUG, 31/05/2010
    else
        B = [];
        dB = [];
        return
    end
else
    % Default, right B-spline
    tt = knots(jmin:jmax+k);
    if issorted(tt)
        [~, col] = histc(x,tt); %#ok
    else
        B = [];
        dB = [];
        return
    end
end

inside = col>=1 & col<=size(M,2);
i = col; % interval index, same size as x
row = find(inside);
col = i(inside); % also i(row)

if k >= 2 && m % i must not be empty, required by GetBreaksArray

    % Construct the array of break indices of sorted x array
    breaks = GetBreaksArray(i, jmax+k);

    %  eqt (5), Carl de Boor "On Calculating with B-spline" 1972 paper
    dt = knots(col+1)-knots(col);
    idt = 1 ./ dt;
    idt(dt == 0) = 0;
    M(row+(col-1)*m) = idt;

    % dM has three dimensions:
        % - first -> abscissa (x), m
        % - second -> basis (j), n
        % - third -> knot derivative, p
    %dM = zeros([size(M) p],cls);
    % Derivative k==0
    ddt = dknots(1,col+1,:)-dknots(1,col,:);   % (1 x m x p)
    ddt = reshape(ddt, [m,1,p]);               % (m x 1 x p)
    ddt(dt == 0,:,:) = 0;
    dM = zeros([size(M)],cls);
    dM(row+(col-1)*m) = idt.*idt;               % (m x n) 
    dM = -ddt.*dM;                              % (m x n x p)

    %%
    for kk=2:k
        % sub-interval to be processed
        jv = jmin:jmax+k-kk;
        % support lengths
        dtv = knots(jv+kk)-knots(jv);
        % recursion loop on sub-intervals
        for c=1:size(jv,2)
            dt = dtv(c);
            if dt~=0
                jj = jv(c);
                % which points x fall inside the support
                r = breaks(jj):breaks(jj+kk)-1; % length mr
                dkleft = dknots(1,jj,:);
                dkright = dknots(1,jj+kk,:);
                ddt = dkright-dkleft;           % ( 1 x 1 x p)
                w = (x(r)-knots(jj)) / dt;      % (mr x 1 x 1)
                dw = -(dkleft + ddt.*w) / dt;   % (mr x 1 x p)
                % eqt (8), Carl de Boor "On Calculating with B-spline" 1972 paper
                Mij = M(r,c);
                Mijp1 = M(r,c+1);
                M(r,c) = w.*Mij + (1-w).*Mijp1;
                dM(r,c,:) = dw.*(Mij - Mijp1) + (w.*dM(r,c,:) + (1-w).*dM(r,c+1,:));
            end
        end
    end

    % Normalized B-Spline basis,
    % B here is noted by N in Carl de Boor "On Calculating with B-spline" 1972 paper
    % definition (4)
    % Loop on sub-intervals
    B = zeros(size(M), cls);
    dB = zeros(size(dM), cls);
    jv = jmin:jmax;
    dtv = knots(jv+k)-knots(jv);
    for c=1:size(jv,2)
        dt = dtv(c);
        if dt~=0
            jj = jv(c);
            % which points x fall inside the support
            r = breaks(jj):breaks(jj+k)-1;
            B(r,c) = M(r,c)*dt;
            ddt = dknots(1,jj+k,:)-dknots(1,jj,:);   % (1 x 1 x p)
            dB(r,c,:) = dM(r,c,:)*dt + M(r,c).*ddt;
        end
    end
else
    B = zeros(size(M), cls);
    B(row+(col-1)*m) = 1;
    dB = zeros([size(B) p], cls);
end

%%
% Map to original vector j
%[tf loc] = ismemberc(j, jmin:jmax); %#ok
if length(j) ~= size(B,2)
    loc = ismembc2(j, jmin:jmax);
    B = B(:,loc);
    dB = dB(:,loc,:);
end

%% Restore the order of original x
if ~isequal(isx, (1:size(B,1)).')
    B = B(isx,:,:);
    dB = dB(isx,:,:);
end

% Basis
n = numel(j);
B = reshape(B, [szx n]);
dB = reshape(dB, [szx n p]);

% Multiply with coefficients to get the spline function
if coefin
    alpha = alpha(:);
    
    B = reshape(B,[],n);
    % Compute function from the coefficients
    B = multMat(B, alpha); % Bug fix 10-Jun-2010
    B = reshape(B, [szx size(alpha,2)]);

    if nargin>=7 % && nargout>=3 % lambda is provided
        % left product with the dual
        lambda = lambda(:).'; % row-vector
        dalpha = lambda*reshape(dB, [szx, n*p]);
        dalpha = reshape(dalpha, [n p]);
    end
    
    % permute the coefficient dimension to the end
    d = length(szx)+1;
    ip = 1:max(ndims(dB),d);
    ip([d end]) = ip([end d]);
    dB = reshape(permute(dB, ip), [], length(alpha));
    dB = multMat(dB, alpha);
    dB = reshape(dB, [szx p]);    

end

end
