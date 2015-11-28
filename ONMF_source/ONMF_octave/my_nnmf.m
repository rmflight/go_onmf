function [wbest,hbest,normbest] = my_nnmf(a,k,varargin)
%NNMF Non-negative matrix factorization.
%   [W,H] = NNMF(A,K) factors the non-negative N-by-M matrix A into
%   non-negative factors W (N-by-K) and H (K-by-M).  The result is not an
%   exact factorization, but W*H is a lower-rank approximation to the
%   original matrix A.  The W and H matrices are chosen to minimize the
%   objective function that is defined as the root mean squared residual
%   between A and the approximation W*H.  This is equivalent to
%
%          D = sqrt(norm(A-W*H,'fro')/(N*M))
%
%   The factorization uses an iterative method starting with random initial
%   values for W and H.  Because the objective function often has local
%   minima, repeated factorizations may yield different W and H values.
%   Sometimes the algorithm converges to solutions of lower rank than K,
%   and this is often an indication that the result is not optimal.
%
%   [W,H,D] = NNMF(...) also returns D, the root mean square residual.
%
%   [W,H] = NNMF(A,K,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following parameter name/value pairs:
% 
%      Parameter    Value
%      'algorithm'  Either 'als' (default) to use an alternating least
%                   squares algorithm, or 'mult' to use a multiplicative
%                   update algorithm.
%      'w0'         An N-by-K matrix to be used as the initial value for W.
%      'h0'         A K-by-M matrix to be used as the initial value for H.
%      'replicates' The number of times to repeat the factorization, using
%                   new random starting values for W and H, except at the
%                   first replication if w0 and h0 are given (default 1).
%                   This tends to be most beneficial with the 'mult'
%                   algorithm.
%      'options'    An options structure as created by the STATSET
%                   function.  NNMF uses these options:  Display, TolX,
%                   TolFun, MaxIter.  Unlike in optimization settings,
%                   reaching MaxIter iterations is treated as convergence.
%
%    Examples:
%        % Non-negative rank-2 approximation of the Fisher iris measurements
%        load fisheriris
%        [w,h] = nnmf(meas,2);
%        gscatter(w(:,1),w(:,2),species);
%        hold on; biplot(max(w(:))*h','VarLabels',{'sl' 'sw' 'pl' 'pw'},'positive',true); hold off;
%        axis([0 12 0 12]);
%
%        % Try a few iterations at several replicates using the
%        % multiplicative algorithm, then continue with more iterations
%        % from the best of these results using alternating least squares
%        x = rand(100,20)*rand(20,50);
%        opt = statset('maxiter',5,'display','final');
%        [w,h] = nnmf(x,5,'rep',10,'opt',opt,'alg','mult');
%        opt = statset('maxiter',1000,'display','final');
%        [w,h] = nnmf(x,5,'w0',w,'h0',h,'opt',opt,'alg','als');
%
%    See also BIPLOT, PRINCOMP, STATSET.

%   Copyright 2007-2009 The MathWorks, Inc. 
%   $Revision: 1.1.6.5 $  $Date: 2009/05/07 18:31:53 $

% Reference:
%   M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
%     Nonnegative Matrix Factorization," Computational Statistics and Data
%     Analysis, vol. 52, no. 1.

% Results are not unique.  This function normalizes W and H so that the
% rows of H have unit length, and the columns of W are ordered by
% decreasing length.

% Check required arguments
[n,m] = size(a);
if ~isscalar(k) || ~isnumeric(k) || k<1 || k>min(m,n) || k~=round(k)
    error('stats:nnmf:BadK',...
          'K must be a positive integer no larger than the number of rows or columns in A.');
end

% Process optional arguments
pnames = {'algorithm' 'w0' 'h0' 'replicates' 'options'};
dflts  = {'als'       []   []   1            []       };
[eid,emsg,alg,w0,h0,tries,options] = ...
        getargs(pnames,dflts,varargin{:});
if ~isempty(eid)
    error(sprintf('stats:nnmf:%s',eid),emsg)
end
    
% Check optional arguments
alg = statgetkeyword(alg,{'mult' 'als'},false,'ALGORITHM','stats:nnmf:BadAlg');
ismult = strncmp('mult',alg,numel(alg));
checkmatrices(a,w0,h0,k);
if ~isscalar(tries) || ~isnumeric(tries) || tries<1 || tries~=round(tries)
    error('stats:nnmf:BadReplicates',...
          'REPLCATES must be a positive integer.');
end

defaultopt = statset('nnmf');
tolx = statget(options,'TolX',defaultopt,'fast');
tolfun = statget(options,'TolFun',defaultopt,'fast');
maxiter = statget(options,'MaxIter',defaultopt,'fast');
dispopt = statget(options,'Display',defaultopt,'fast');
dispnum = strmatch(lower(dispopt), {'off','notify','final','iter'}) - 1;

% Special case, if K is full rank we know the answer
if isempty(w0) && isempty(h0)
    if k==m
        w0 = a;
        h0 = eye(k);
    elseif k==n
        w0 = eye(k);
        h0 = a;
    end
end

% Prepare to start iterations
normbest = Inf;
ws = warning('off','MATLAB:rankDeficientMatrix');
if dispnum > 1 % 'iter' or 'final'
    fprintf('    rep\titeration\t   rms resid\t   |delta x|\n');
end

try
    for j=1:tries
        % Get random starting values if required
        if isempty(w0) || j>1
            w0 = rand(n,k);
        end
        if isempty(h0) || j>1
            h0 = rand(k,m);
        end

        % Perform a factorization
        [w,h,dnorm,iters] = nnmf1(a,w0,h0,ismult,maxiter,tolfun,tolx,...
                                  dispnum,j);

        % Save if this is the best so far
        if dnorm<normbest
            wbest = w;
            hbest = h;
            normbest = dnorm;
        end
    end
    if dispnum > 1   % 'final' or 'iter'
       fprintf('Final root mean square residual = %g\n',normbest);
    end
catch ME
    warning(ws);
    rethrow(ME);
end
warning(ws);

if normbest==Inf
    error('stats:nnmf:NoSolution',...
          'Algorithm could not converge to a finite solution.')
end

% Put the outputs in a standard form - first normalize h
hlen = sqrt(sum(hbest.^2,2));
if any(hlen==0)
    warning('stats:nnmf:LowRank',...
            'Algorithm converged to a solution of rank %d rather than %d as specified.',...
            k-sum(hlen==0), k);
    hlen(hlen==0) = 1;
end
wbest = bsxfun(@times,wbest,hlen');
hbest = bsxfun(@times,hbest,1./hlen);

% Then order by w
[ignore,idx] = sort(sum(wbest.^2,1),'descend');
wbest = wbest(:,idx);
hbest = hbest(idx,:);
end

% -------------------
function [w,h,dnorm,iters] = nnmf1(a,w0,h0,ismult,maxiter,tolfun,tolx,...
                                   dispnum,repnum)
% Single non-negative matrix factorization
nm = numel(a);
sqrteps = sqrt(eps);
iters = maxiter;
dispfmt = '%7d\t%8d\t%12g\t%12g\n';

for j=1:maxiter
    if ismult
        % Multiplicative update formula
        numer = w0'*a;
        h = h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer)));
        numer = a*h';
        w = w0 .* (numer ./ (w0*(h*h') + eps(numer)));
        % w = w0 .* (numer ./ ((w0*h*a'*w0) + eps(numer)));
        % w = w0 .* (numer ./ ((mtimesx(w0, mtimesx(h, mtimesx(a',w0,'SPEEDOMP'),'SPEEDOMP'),'SPEEDOMP')) + eps(numer)));
    	% diag_w = diag(sum(w));
        % w = mtimesx(w, inv(diag_w), 'SPEEDOMP');
    else
        % Alternating least squares
        h = max(0, w0\a);
        w = max(0, a/h);
    end
    
    % Get norm of difference and max change in factors
    d = a - w*h;
    dnorm = sqrt(sum(sum(d.^2))/nm);
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
    
    % Check for convergence
    if j>1
        if delta <= tolx
            iters = j;
            break;
        elseif dnorm <= tolfun*dnorm0
            iters = j;
            break;
        elseif j==maxiter
            break
        end
    end

    if dispnum>2 % 'iter'
        fprintf(dispfmt,repnum,j,dnorm,delta);
    end

    % Remember previous iteration results
    dnorm0 = dnorm;
    w0 = w;
    h0 = h;
end

if dispnum>1   % 'final' or 'iter'
    fprintf(dispfmt,repnum,j,dnorm,delta);
end

end

% ---------------------------
function checkmatrices(a,w,h,k)
% check for non-negative matrices of the proper size

if ndims(a)~=2 || ~isnumeric(a) || ~isreal(a) || any(any(a<0)) || any(any(~isfinite(a)))
    error('stats:nnmf:BadA','A must be a matrix of non-negative values.')
end
[n,m] = size(a);
if ~isempty(w)
    if ndims(w)~=2 || ~isnumeric(w)|| ~isreal(w) || any(any(w<0)) || any(any(~isfinite(w)))
        error('stats:nnmf:BadW','W must be a matrix of non-negative values.')
    elseif ~isequal(size(w),[n k])
        error('stats:nnmf:BadW',...
              'The size of the W matrix must be %d-by-%d.',n,k)
    end
end
if ~isempty(h)
    if ndims(h)~=2 || ~isnumeric(h)|| ~isreal(h) || any(any(h<0)) || any(any(~isfinite(h)))
        error('stats:nnmf:BadH','H must be a matrix of non-negative values.')
    elseif ~isequal(size(h),[k m])
        error('stats:nnmf:BadH',...
              'The size of the H matrix must be %d-by-%d.',k,m)
    end
end
end % function


