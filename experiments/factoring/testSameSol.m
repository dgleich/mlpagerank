% This script aims to test equivalence of the solutions to tensor 
% PageRank model, 2-step tensor PageRank model, and higher order 
% PageRank models under the condition that R = kron(e', Q);

%% Load data
load('sameSol.mat');
n = size(R,1);
%% Compute PageRank
P = R(1:n,1:n); xpr = (eye(n) - alpha*P)\v; xpr = xpr/sum(xpr);
%% Compute tensor PageRank
x1 = tensorRank(alpha, R, v,v); % x1 -- solution to tensor PageRank model
%% Compare
[xpr x1 xpr-x1]
%% Compute the higher order solutions
[X1, X2] = oneTwoStepChain(alpha, R, v); % X1 -- solution to 1-step higher order PageRank model
                                         % X2 -- solution to 2-step higher order PageRank model
%% Compare
fprintf('X1\n');
[X1 kron(x1,x1) X1-kron(x1,x1)]
fprintf('x2\n');
[X1 kron(x1,x1) X1-kron(x1,x1)]
%% X1 = X2, but they are not kronecker products...
X1 = reshape(X1,n,n);
x21 = sum(X1,2)
[x1 x21 x21-x1]
%% So the marginals' hold.
