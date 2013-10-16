function J = getJ(alpha, R, v, x)
% J = getJ(alpha, R, v, x)
% function getJ computes the Jacobian of
% [I - alpha/2*R*(kron(x, I) + kron(I, x))]^(-1)*(1-alpha)*v = x
% with respect to x
n = size(R, 1); % dimension
J = zeros(n, n);
I = eye(n);
% compute the deriverative to each component x_i
A = I - alpha/2*R*(kron(x, I) + kron(I, x));
A1 = inv(A);
% compute A^-1
% A1 = zeros(n,n);
% for j = 1:n
%     e = zeros(n,1);
%     e(j) = 1;
%     A1(:,j) = A \ e; 
% end
for i = 1:n
    B = zeros(n^2, n);
    B((i-1)*n + 1:i*n, :) = I;
    C = zeros(n^2, n);
    for k = 1:n
        C((k-1)*n + i, k) = 1;
    end
    J(:,i) = (1-alpha)*alpha/2*A1*R*(B+C)*A1*v;
end
J = J - I;
end