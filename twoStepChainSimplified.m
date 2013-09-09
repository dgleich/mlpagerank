function [xtrue, xhist] = twoStepChainSimplified(alpha, R, v, x)
% [xtrue, xhist] = twoStepChainSimplified(alpha, R, v, x)
% The function is simplified by replacing X with kron(x,x). The goal is to 
% compare the solution with those from oneTwoStepChain.m.
% The equaiton of the two-step Markov chain becomes,
% i.e., x =  [alpha^2*R*P + alpha*(1-alpha)*R*V]*kron(x,x) + (1-alpha)*v
 
P = convertR2P(R);
I = eye(size(v,1));
e = ones(size(v));
V = kron(e',kron(I, v));
tol = 1e-12;
niter = 10000;
xhist = zeros(size(x,1), niter);
xcur = x;
for n = 1:niter
    M = alpha*alpha*R*P + alpha*(1-alpha)*R*V;
    xn = M*kron(xcur,xcur) + (1-alpha)*v;
    xn = xn ./ norm(xn, 1);
    
    xhist(:, n) = xn;
    if norm(xn - xcur, inf) <= tol
        xtrue = xn;
        xhist = xhist(:, 1:n);
        break;
    else
        xcur = xn;
    end
end
if size(xhist, 2) == niter
    fprintf('reaching max iter times!\n');
    xtrue = xhist(:, niter);
end
% test if the solution is the true one
r = (alpha*alpha*R*P + alpha*(1-alpha)*R*V)*kron(xtrue, xtrue) + (1-alpha)*v - xtrue;
if norm(r, 1) <= tol
    fprintf('True solution\n');
end
end