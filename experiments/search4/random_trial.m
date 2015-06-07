function [Rf,niters,nnzs] = random_trial(R,alpha)

n = size(R,1);
v = ones(n,1)/n;

startnz = sum(R(:));

nnzs = zeros(startnz,1);
niters = zeros(startnz,1);

Rf = {};

for i=1:startnz
    % compute tpr
    [x,hist,flag] = tpr4_power(normout(R')', ...
                        alpha, v, 'maxiter', 1e5);
                    
    if ~flag
        Rf = {R};
        nnzs = nnzs(1:i-1);
        niters = niters(1:i-1);
        return
    else
        nnzs(i) = sum(R(:));
        niters(i) = length(hist);
         
        % find the non-zero with the largest value 
        % and remove it
        X = reshape(kron(kron(x,x),kron(x,x)), n, n, n, n);
        T = reshape(R, n, n, n, n);
        V = X.*(1./T);
        [~,ind] = min(V(:));
        assert(R(ind) > 0)
        R(ind) = 0;
    end
end

assert( any(R(:)) == 0)
