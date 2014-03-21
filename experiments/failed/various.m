%%
load test_case/6x6x6_not_converge_non_shift.mat
xbad = x1;
alpha = 0.99;

n = size(R1,1);
x2 = oneTwoStepChain(alpha,R1,ones(6,1)/6);
X = reshape(x2,n,n);
Sx = sum(X,2);

%%
P = tensorpr3(R1,alpha);
%%
P.residual(Sx)
J = P.jacobian(Sx,0.5);
eig(J)

%%
P.residual(x1)
J = P.jacobian(x1,1);
eig(J)

%% Okay, David make a mistake in thinking of these things. Oops! But it led
% to some useful code.

%% New ideas
P = sparse(convertR2P(R1));
spy(P)
cc = components(P)
test_dag(P)
%%
[p,p,b,b] = dmperm(P + speye(size(R1,2)));
spy(P(p,p)+speye(size(R1,2)))

