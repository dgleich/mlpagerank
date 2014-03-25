function P = tensor2markov(R)
% TENSOR2MARKOV Convert the transition tensor to a 2nd order Markov matrix
%
% P = tensor2markov(R) returns the n^2 by n^2 transition matrix for the
% second order Markov chain encoding the dynamics:
%    prob( X_t+1 = i | X_t = j, X_t-1 = k ) = P(i,j,k)
% where i,j,k vary between 1 to n.

n = size(R, 1);
P = zeros(n^2, n^2);
for i = 1:n     % group i
    tmp = zeros(n^2, n);
    for j = 1:n % column j
        ej = zeros(n, 1);
        ej(j) = 1;
        tmp(:, j) = kron(ej, R(:, (i-1)*n + j));
    end
    P(:, (i-1)*n +1: i*n) = tmp;
end
end