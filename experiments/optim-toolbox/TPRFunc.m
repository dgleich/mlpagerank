function [f,J] = TPRFunc(tpr,x)

f = tpr.residual(x);
J = tpr.jacobian(x) - eye(size(x,1));
%FuncJacobian = @(x) deal(Func(x), Jacobian(x));
