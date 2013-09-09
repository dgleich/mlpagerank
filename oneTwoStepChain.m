function [X1, X2] = oneTwoStepChain(alpha, R, v)
% [X1, X2] = oneTwoStepChain(alpha, R, v)
% this function computes the true eigenvalues of the one-step Markov chain
% and the two-step Markov chain for the following equations, i.e.,
% [alpha*P + (1-alpha)*V]*X = X  and 
% [alpha*P + (1-alpha)*V]*[alpha*P + (1-alpha)*V]*X = X
% where P = convertR2P(R) and V = kron(v, kron(I, e')).
tol = 1e-12;
v = v ./ sum(v);
P = convertR2P(R);
e = ones(size(v));
I = eye(size(v,1));
V = kron(e', kron(I,v));
Q = alpha*P + (1-alpha)*V;
% X1 is the eigen vector of Q with eigenvalue = 1
X1 = null(Q - eye(size(Q)));
X1 = X1 / sum(X1);
% test validation
if norm(X1 - Q*X1, 1) <= tol
    fprintf('X1 is the true solution for one-step chain\n');
end
Q = mpower(Q, 2);
% X2 is the eigen vector of Q^2 with eigenvalue = 1
X2 = null(Q - eye(size(Q)));
X2 = X2 ./ sum(X2);
if norm(X2 - Q*X2, 1) <= tol
    fprintf('X2 is the true solution for two-step chain\n');
end
end