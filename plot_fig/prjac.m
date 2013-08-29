function J = prjac(x,alpha,R,v,gamma)

n = size(R,1);
J = alpha*gamma*R*(kron(x,eye(3)) + kron(eye(3),x)) + (1-gamma)*eye(3);

