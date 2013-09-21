% This script aims to test equivalence of the solutions to tensor PageRank model,
% 2-step tensor PageRank model, and higer order PageRank models
% under the condition that R = kron(e', Q);

load('sameSol.mat');
x1 = tensorRank(alpha, R, v,v); % x1 -- solution to tensor PageRank model
x2 = twoStepChainSimplified(alpha, R, v,v); % x2 -- solution to 2-step tensor PageRank model
[X1, X2] = oneTwoStepChain(alpha, R, v); % X1 -- solution to 1-step higher order PageRank model
                                         % X2 -- solution to 2-step higher order PageRank model
fprintf('norm(x1-x2,1) is %f\n', norm(x1-x2, 1));
fprintf('norm(X1-X2,1) is %f\n', norm(X1-X2, 1));
% test if X1 == kron(x1, x1)
fprintf('norm(X1-kron(x1,x1), 1) is %f\n', norm(X1-kron(x1,x1), 1));