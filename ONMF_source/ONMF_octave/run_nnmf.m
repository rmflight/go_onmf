function [w,h,dnorm,iters] = run_nnmf(a,w0,h0,ismult,maxiter,tolfun,tolx,...
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
        h = h0 .* (numer ./ (((w0'*w0)*h0) + eps(numer)));

        numer = a * h';
        % w = w0 .* (numer ./ (w0*(h*h') + eps(numer)));
        w = w0 .* (numer ./ ((((w0*h)*a')*w0) + eps(numer)));
    	diag_w = diag(sum(w));
        % w = mtimesx(w, inv(diag_w), 'SPEEDOMP');
	w = w * inv(diag_w);
    else
        % Alternating least squares
        h = max(0, w0\a);
        w = max(0, a/h);
    end
    
    % Get norm of difference and max change in factors
    d = a - w * h;
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

    % Remember previous iteration results
    dnorm0 = dnorm;
    w0 = w;
    h0 = h;
end
