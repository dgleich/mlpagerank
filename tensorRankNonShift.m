function [xtrue, xhist, kappa, flag] = tensorRankNonShift(alpha, R, v, x)
% [xtrue, xhist, kappa, flag] = tensorRankNonShift(alpha, R, v, x)
% new update for the tensorRank problem
% x^(k+1) = (1-alpha)*(I - alpha/2*R*(kron(x^(k), I) + kron(I,
% x^(k))))^(-1)*v
% kappa records the ratio of ||x^(k+1) - x^(k)||_1 / ||x^(k) - x^(k-1)||_1
tol = 1e-12;
niter = 10000;
xhist = zeros(size(x,1), niter);
kappa = zeros(niter - 2, 1);
xcur = x;
I = eye(max(size(x)));
for n = 1:niter
    A = kron(xcur, I) + kron(I, xcur);
    A = I - alpha/2 .*R*A;
    b = (1-alpha).*v;
    xn = A \ b;
    xn = xn ./ norm(xn, 1);
    
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
flag = valSol(alpha, R, xtrue, v);
end
