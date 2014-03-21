%% Run the power method of Li and Ng on a PageRank modified tensor
load('./test_case/ling1') % load the first example for li and ng
tpr3 = tensorpr3(P, 1);
[~,R] = tpr3.markov();
% li-ng's method
x0 = ones(size(R,1), 1)/size(R,1);
xt = x0;
for i=1:100
    x = R*kron(xt,xt);
    if norm(x - xt, 1) <= 1e-8
        break;
    else
        xt = x;
    end
end
y = tpr3.solve();
if norm(x-y,1) <= 1e-8
    fprintf('ling1 test successful.\n');
else
    fprintf('ling1 test failed.\n');
end