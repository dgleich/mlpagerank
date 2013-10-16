function test_tensorpr3
%%
R = kron(ones(1,3),eye(3));
x = tensorpr3(R,0.85).solve();
assert(norm(x - 1/3)<10*eps(1),'failed test on uniform');

%%
R = kron(ones(1,3),eye(3));
x = tensorpr3(R,0.85,[1/2;1/3;1/6]).solve();
assert(norm(x - [1/2;1/3;1/6])<10*eps(1),'failed test on spec');
%%
P = rand(3,3);
P = P*diag(1./sum(P));
R = kron(ones(1,3),P);
xpr = (eye(3) - 0.85*P)\((0.15)*1/3*ones(3,1));
TPR = tensorpr3(R,0.85);
x = TPR.solve();
assert(norm(x - xpr)<1e-8, 'failed on prtest');
x = TPR.solve('tol',eps(1));
assert(norm(x - xpr)<10*eps(1), 'failed on high-precision prtest');