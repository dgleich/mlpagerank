function flag =  valSol(alpha, R, x, v)
% flag =  valSol(alpha, R, x, v)
% this function test if x = alpha * R * kron(x, x) + (1-alpha)*v
res = alpha .* R * kron(x, x) + (1-alpha) .* v;
a = res - x;
if norm(a, 1) <= 1e-10
    fprintf('True solution.\n');
    flag = 0;
else
    fprintf('False solution!\n');
    flag = 1;
end
%fprintf('The 1-norm is ');
%disp(norm(a,1));
end
