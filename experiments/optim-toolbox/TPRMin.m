function [f,g] = TPRMin(tpr,x)

% we want to minimize
% Let F(x) be the vector residual
% we want to minimize
% f(x) = 1/2*F(x)'*F(x)
% so the gradient is
% g(x) = 2*J(x)

F = tpr.residual(x);
f = 0.5*F'*F;
g = (tpr.jacobian(x) - eye(size(x,1)))'*F;

