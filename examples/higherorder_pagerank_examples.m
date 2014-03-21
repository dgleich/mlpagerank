%% Example for paper

load examples/example1
n = size(R,1);

v = ones(3,1)/n;
alpha = 0.85;

V = kron(kron(ones(1,n),eye(n)),v);
[X,D] = eig((alpha*P + (1-alpha)*V));

x = X(:,1);
x = x/sum(x);

X = reshape(x,3,3)

