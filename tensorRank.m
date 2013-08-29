function [xtrue, xhist, flag, sumx, diffsumx, kappa, flagk] ...
    = tensorRank(alpha, R, v, x, gamma)
% [xtrue, xhist, flag, sumx, diffsumx, kappa, flagr] = tensorRank(alpha, R, v, x, gamma)
% xtrue is the true value of vector x (computed after nth iteration)
% xhist records each iteration
% sumx records the sum of vector x in each iteration
% diffsumx records the difference between x_n and 1
% kappa records the ratio of ||x^(k+1) - x^(k)||_1 / ||x^(k) - x^(k-1)||_1
% flag -- 0 converge;  flag -- 1, not converge
% gamma -- shifting factor, default value is 0.5
% flagk -- indicator for kappa > 1
if nargin < 5
	gamma = 0.5;
end
tol = 1e-12;
niter = 10000;
flag = 0; 
flagk = 0;
xhist = zeros(size(x,1), niter);
sumx = zeros(niter, 1);
diffsumx = zeros(niter, 1);
kappa = zeros(niter - 2, 1);
xcur = x;
Gamma = 1 / (1+gamma);
for n = 1 : niter
    %xn = (alpha*R*kron(xcur, xcur) + ...
    %    (1-alpha)*v + gamma*xcur)./(1+gamma);
    %xn = xn ./ norm(xn, 1); % normalization step
    y = alpha*R*kron(xcur, xcur);
    z = y * Gamma + Gamma*(1-sum(y))*v;
    xn = z + (1-sum(z))*xcur;
    sumx(n) = sum(xn);
    diffsumx(n) = norm(1-sum(xn), 1);
    
    xhist(:, n) = xn;
    if n > 2
            kappa(n - 2) = norm(xhist(:,n) - xhist(:,n-1), 1) / norm(xhist(:,n-1) - xhist(:,n-2), 1);
			if kappa(n - 2) > 1
				flagk = 1;
			end
    end
    if norm(xn - xcur, inf) <= tol
        xtrue = xn;
        xhist = xhist(:, 1:n);
        sumx = sumx(1 : n);
        diffsumx = diffsumx(1:n);
        kappa = kappa(1:n-2);
        break;
    else
        xcur = xn;        
    end
end
if size(xhist, 2) == niter
    fprintf('Reaching the max iteration times, not converge!\n');
    xtrue = xhist(:, niter);
    flag = 1;
end
valSol(alpha, R, xtrue, v);
end
