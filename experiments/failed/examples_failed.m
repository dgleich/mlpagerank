%% Load data
load('test_case/6x6x6_not_converge_non_shift');

%%
R = R1;
xtrue = x1;
alpha = 0.99;
v = ones(6,1)/6;
[x,xhist] = tensorRank(alpha,R,v,v);

%%
plot(xhist');

%% 
Z = xhist - xtrue*ones(1,length(xhist)); 
Z2 = sum(Z.*Z);
plot(Z2(1:200));

%%
R = R3;
alpha = 0.99;
v = ones(6,1)/6;
[x,xhist] = tensorRank(alpha,R,v,v);

%% 
% But it works with more iterations
R = R3;
alpha = 0.99;
v = ones(6,1)/6;
[x,xhist] = tensorRank(alpha,R,v,v,0.5,1e5);
