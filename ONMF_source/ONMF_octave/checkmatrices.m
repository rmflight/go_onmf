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