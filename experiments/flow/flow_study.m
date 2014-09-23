addpath('~/dev/matlab-bgl');

%%


n = 2;
search_tensors(n,@dynsys_converge,sprintf('flow-%i',n),'filter',@strong_connectivity_filter);

%%
load flow-2-found

T = T1;
n = size(T,1);
alpha = 1;
v = ones(n,1)/n;
R = reshape(T,n,n^2);
R = normout(R')';
tpr = tensorpr3(R,alpha,v);
[xhist,thist,rhist] = tpr.dynsys('method','ode45');

%%
plot(xhist')