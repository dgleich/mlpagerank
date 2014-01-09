%% Run the power method of Li and Ng on a PageRank modified tensor
load('./test_case/ling1') % load the first example for li and ng
tpr3 = tensorpr3(P, 1);
[~,R] = tpr3.markov();
x0 = ones(size(R,1), 1)/size(R,1);
xt = x0;
for i=1:100
    x = R*kron(xt,xt);
    if norm(x - xt, 1) < 1e-8
        break;
    else
        xt = x;
    end
end