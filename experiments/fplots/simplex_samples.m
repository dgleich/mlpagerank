function [T,Z,X,C]=simplex_samples(n)
% 
% [T,S] = simplex_samples(n)
%   S - points on the simplex
%   T - points on the triangle (for visualization)
%

ts = 2*pi.*(1:3)/3 - pi/6;
X = [cos(ts); sin(ts)];

T = triangle_grid(n,X);

% Compute simplex coordinates for each point in T.
C = [X(:,1) - X(:,3)  X(:,2) - X(:,3)];
Z = C\(T - X(:,3)*ones(1,size(T,2)));
Z(abs(Z) < 2*eps(1)) = 0;
Z = [Z; 1-sum(Z)];
Z(abs(Z) < 2*eps(1)) = 0;
for i=1:size(Z,2)
    Z(:,i) = Z(:,i)/sum(Z(:,i));
end


