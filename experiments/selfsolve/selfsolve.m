% Did we ever try the following?
% 
% If we solve 
% 
% alpha R x^2 + (1-alpha) v = x
% 
% and then use the solution x as the v for the next iterate.
% 
% If it converges, then it must converge to:
% alpha R x^2 + (1-alpha) x = x
% or
% alpha R x^2 + x - alpha x = x
% or
% alpha R x^2 - alpha x = 0
% or
% R x^2 = x
% which means if we set R to be a PageRank system itself,
% then we MIGHT get convergence here 

%% Load a bad example
load('../../tensors/4x4x4_not_converge_non_shift');

%%
R = R16;
n = size(R,1);
alpha = 0.99;
v = ones(n,1)./n;
T = tensorpr3(R,alpha,v);
x = T.solve();

%%
Rt = alpha*R + (1-alpha)*v*ones(1,n^2);
at = alpha/2;
xt = v;
for i=1:1000
    Tt = tensorpr3(Rt,at,xt);
    xt2 = Tt.solve();
    fptrintf('%.4e\n', norm(xt-xt2,1));
    xt = xt2;
end
%%
T.residual(xt)