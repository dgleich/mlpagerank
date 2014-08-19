R = [0 0 0 0 0 0 1/3 1 0;
     0 0 0 0 1 0 1/3 0 1;
     1 1 1 1 0 1 1/3 0 0];
v = [0; 1; 0];
alpha = 0.99;

tpr = tensorpr3(R,alpha,v);
x = tpr.shifted()
y = tpr.shifted('randinit', true)

%%
M = tpr.markov;

%% Try a perturbation
v1 = v + 1e-8*rand(3,1);
v1 = v1/sum(v1);
tpr = tensorpr3(R,alpha,v1);
tpr.residual(v)

%%
x = tpr.shifted()
y = tpr.shifted('randinit', true)

