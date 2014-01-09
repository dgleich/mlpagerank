%% Run the power method of Li and Ng on a PageRank modified tensor
load ling1 % load the first example for li and ng
[~,R] = tensorpr3(P,1).markov();
x0 = ones(size(R,1))/size(R,1);
x = x0;
for i=1:100
    x = R*kron(x,x);
end