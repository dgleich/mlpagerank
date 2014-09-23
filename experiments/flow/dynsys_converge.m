function rval = dynsys_converge(T)
n = size(T,1);
alpha = 1;
v = ones(n,1)/n;
R = reshape(T,n,n^2);
R = normout(R')';
tpr = tensorpr3(R,alpha,v);
[xhist,thist,rhist] = tpr.dynsys('method','ode45');
if rhist(end) > 1e-6
    rval = 1;
else
    rval = 0;
end