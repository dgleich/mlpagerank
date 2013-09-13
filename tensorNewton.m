function [xtrue, xhist, kappa] = tensorNewton(alpha, R, v, x, maxIter, ntol)
% function [xtrue, xhist, kappa] = tensorNewton(alpha, R, v, x, maxIter, ntol)
% Newton's method for the tensor problem 
% [I - alpha*R*(kron(x_n, I) + kron(I, x_n))]*x_(n+1) = ...
%  -alpha*R*kron(x_n,x_n) + (1-alpha)*v
% kappa records the ratio of ||x^(k+1) - x^(k)||_1 / ||x^(k) - x^(k-1)||_1
if nargin < 5
    niter = 1e4;
    tol = 1e-12;
else
    niter = maxIter;
    if nargin < 6
        niter = maxIter;
    else
        tol = ntol;
    end
end

xhist = zeros(size(x,1), niter);
kappa = zeros(niter - 2, 1);
xcur = x;
I = eye(max(size(x)));
for n = 1:niter
    A = I - alpha.*R*(kron(xcur, I) + kron(I, xcur));
    b = -alpha.*R*kron(xcur, xcur) + (1-alpha)*v;
    xn = A \ b;
    xn = xn ./ sum(xn);
    
    xhist(:, n) = xn;
    if n > 2
        kappa(n - 2) = norm(xhist(:,n) - xhist(:,n-1), 1) / norm(xhist(:,n-1) - xhist(:,n-2), 1);
    end
    if norm(xn - xcur, inf) <= tol
        xtrue = xn;
        xhist = xhist(:, 1:n);
        kappa = kappa(1:n-2);
        break;
    else
        xcur = xn;
    end
end
if size(xhist, 2) == niter
    fprintf('reaching max iter times!\n');
    xtrue = xhist(:, niter);
end
valSol(alpha, R, xtrue, v);
end