%% Load a bad example
load('../test_case/4x4x4_not_converge_non_shift', 'R10');
%%
load('../test_case/6x6x6_not_converge_non_shift');
%%
R = R1;
n = size(R,1);
alpha = 0.99;
v = ones(n,1)./n;
T = tensorpr3(R,alpha,v);
x = T.solve();

%%
R=R10;
n = size(R,1);
alpha = 0.99;
v = ones(n,1)./n;
Rt = alpha*R + (1-alpha)*v*ones(1,n^2);
at = alpha/2;
xt = v;
xhist = zeros(n, 1000);
for i=1:1000
Tt = tensorpr3(Rt,at,xt);
xt2 = Tt.solve();
xhist(:,i) = xt2;
%fprintf('%.4e\n', norm(xt-xt2,1));
xt = xt2;
end
plot([v xhist]');