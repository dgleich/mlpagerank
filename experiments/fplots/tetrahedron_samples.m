function [T,X,Z,C] = tetrahedron_samples(n)

T = [1 0 0 1; 0 1 0 1; 0 0 1 1];
Z = tetrahedron_grid(n, T, 1); % just let Matlab allocate it
%% vertices at the origin
C = [T(:, 1) - T(:,4) T(:,2)-T(:,4) T(:,3)-T(:,4)];
% map the sample points to the 4x4x4 tensor elements
X = C\(Z - T(:,4)*ones(1, size(Z,2)));
X(abs(X) < 2*eps(1)) = 0;
X = [X; 1-sum(X)];
X(abs(X) < 2*eps(1)) = 0;
for i=1:size(X,2)
    X(:,i) = X(:,i)/sum(X(:,i));
end