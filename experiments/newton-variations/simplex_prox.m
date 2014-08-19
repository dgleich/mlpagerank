function xn = simplex_prox(x)
% SIMPLEX_PROX Compute the proximal point on the probability simplex
%
% Given any vector x, xn = simplex_prox(x) computes the closest point in a
% least-squares sense to x that is within the probability simplex. 
%
% See section 6.2.5 in https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf
%


% The solution xn = max(x - nu,0) such that the resulting vector xn has sum
% 1.

% To find this, we just do bisection on nu.

% The region for nu is max(x)-1, to max(x) because if max(x) is really
% large, then we just project onto the largest coordinate and rescale.
nu = bisection(@(nu) sum(max(x - nu,0)) - 1, max(x) -1, max(x), eps);
xn = max(x - nu,0);

function x0 = bisection(f, a, b, delta)
% BISECTION Find a point where f(x) = 0 through bisection
% x0 = bisection(f, a, b, delta) does an interval bisection search
% to find a region of size delta that contains a zero of
% the function f, by default delta = 2.2e-16, the machine eps.
fa = f(a); fb = f(b); assert(sign(fa*fb) <= 0); maxit = 52;
for i=1:maxit
    ab2 = 0.5*a + 0.5*b; fab2 = f(ab2); if abs(fab2) < eps, break; end
    if abs(b-a) <= delta, break; end
    if sign(fab2*fb) <= 0, a = ab2; fa = fab2; 
    else b = ab2; fb = fab2; end
end
x0 = ab2;