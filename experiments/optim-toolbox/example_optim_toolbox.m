%% Here is an exmaple where the standard fsolve function fails
M = load('../../tensors/mtm4');
alpha = 0.99;
R = M.R4_12;
tpr = tensorpr3(R,alpha);
n = size(R,1);
x0 = tpr.v;

% Can't be fancy because deal won't work with partial outputs :(
%Func = @(x) tpr.residual(x);
%Jacobian = @(x) tpr.jacobian(x) - eye(n);
%FuncJacobian = @(x) deal(Func(x), Jacobian(x));

FuncJacobian = @(x) TPRFunc(tpr,x);

options=optimset('Display','iter','Jacobian','on','DerivativeCheck','on');
[x,F,exitflag,output,JAC] = fsolve(FuncJacobian,x0,options);

sum(x)
x

%%

options=optimset('Display','iter','Jacobian','on','Algorithm','trust-region-reflective');
[x,F,exitflag,output,JAC] = fsolve(FuncJacobian,x0,options);

sum(x)
x

%%
options=optimset('Display','iter','Jacobian','on','Algorithm','levenberg-marquardt');
[x,F,exitflag,output,JAC] = fsolve(FuncJacobian,x0,options);

sum(x)
x

%% So none of Matlab's solvers can work on this case.
% What about treating this as a constrained minimization problem?

M = load('../../tensors/mtm4');
alpha = 0.99;
R = M.R4_17;
xtrue = M.R4_Properties.R4_17.alpha99.sols;
tpr = tensorpr3(R,alpha);
n = size(R,1);
x0 = tpr.v;

fungrad = @(x) TPRMin(tpr,x);

options=optimset('GradObj','on','Display','iter','DerivativeCheck', 'on', 'TolFun',1e-12,'MaxFunEval',10000,'MaxIter',10000);
[x,fval,exitflag,output] = fminunc(fungrad,x0,options);
sum(x)
[x xtrue]

%%
M = load('../../tensors/mtm6');
alpha = 0.99;
R = M.R6_1;
tpr = tensorpr3(R,alpha);
n = size(R,1);
x0 = tpr.v;

fungrad = @(x) TPRMin(tpr,x);

options=optimset('GradObj','on','Display','iter','DerivativeCheck', 'on', 'TolFun',1e-14,'MaxFunEval',10000,'MaxIter',10000);
[x,fval,exitflag,output] = fminunc(fungrad,x0,options);
[sum(x) norm(tpr.residual(x),1)]
[x tpr.residual(x)]

%% This algorithm seems to be fairly effective.
M = load('../../tensors/mtm4');
alpha = 0.9;
R = M.R4_19;
tpr = tensorpr3(R,alpha);
n = size(R,1);
x0 = tpr.v;

fungrad = @(x) TPRMin(tpr,x);

options=optimset('GradObj','on','Display','iter','DerivativeCheck', 'on', 'TolFun',1e-14,'MaxFunEval',10000,'MaxIter',10000);
[x,fval,exitflag,output] = fmincon(fungrad,x0,[],[],ones(1,n),1,0,1,[],options);
[sum(x) norm(tpr.residual(x),1)]
[x tpr.residual(x)]

