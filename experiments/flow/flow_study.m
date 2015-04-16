addpath('~/dev/matlab-bgl');

%%


n = 2;
search_tensors(n,@dynsys_converge,sprintf('flow-%i',n),'filter',@strong_connectivity_filter);

%%
n = 2;
search_tensors(n,@dynsys_converge,sprintf('flow-%i',n),'filter',@strong_connectivity_filter,'max',10);

%% we now have a theorem that shows that n=2 has no counter examples

%% so now let's look at n=3;
n = 3;
search_tensors(n,@dynsys_converge,sprintf('flow-%i',n),'filter',@strong_connectivity_filter);



%%
load flow-3-found

T = T2;
n = size(T,1);
alpha = 1;
v = ones(n,1)/n;
R = reshape(T,n,n^2);
R = normout(R')';
tpr = tensorpr3(R,alpha,v);
[xhist,thist,rhist] = tpr.dynsys('method','ode45','maxiter',100000);

%%
plot(xhist')