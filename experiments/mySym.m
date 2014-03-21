function [Y,Yall] = mySym(R)
syms x1 x2 x3 x4 
alpha = 0.99;
v = 1/4;
[sol_x1, sol_x2, sol_x3, sol_x4] = ...
	solve(x1 == alpha*(0 + R(1,13)*x4*x1 + R(1,16)*x4*x4) + (1-alpha)*v, ...
		  x2 == alpha*(0 + R(2,6)*x2*x2 + R(2,8)*x2*x4 + R(2,10)*x3*x2 + R(2,14)*x4*x2) + (1-alpha)*v, ...
		  x3 == alpha*(0 + R(3,7)*x2*x3 + R(3,10)*x3*x2 + R(3,11)*x3*x3) + (1-alpha)*v, ...
		  x4 == alpha*(0 + R(4,1)*x1*x1 + R(4,2)*x1*x2 + R(4,3)*x1*x3 + R(4,4)*x1*x4 + R(4,5)*x2*x1 + R(4,9)*x3*x1 + R(4,12)*x3*x4 + R(4,13)*x4*x1 + R(4,15)*x4*x3) + (1-alpha)*v);
Y = double([sol_x1 sol_x2 sol_x3 sol_x4 ]');
k = 1;
X = zeros(size(Y));
for j = 1:size(Y, 2)
	if (min(Y(:, j)) >= 0 && sum(Y(:, j)) == 1)
		X(:, k) = Y(:, j);
		k = k + 1;
	end
end
Yall = Y;
Y = X(:, 1:k-1);
end
