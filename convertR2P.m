function P = convertR2P(R)
% P = convertR2P(R)
% the function convert R(nxn^2) to P(n^2 x n^2)
n = size(R, 1);
P = zeros(n^2, n^2);
for i = 1:n     % group i
    tmp = zeros(n^2, n);
    for j = 1:n % column j
        ei = zeros(n, 1);
        ei(i) = 1;
        tmp(:, j) = kron(R(:, (i-1)*n + j), ei);
    end
    P(:, (i-1)*n +1: i*n) = tmp;
end
end