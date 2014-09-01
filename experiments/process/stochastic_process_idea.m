%% Test a wacky stochastic process that may model this vector

load ../../tensors/mtm4
alpha = 0.25; % we want fast convergence if this is the case.

R = R4_1;
%%

n = size(R,1);
P = reshape(R,n,n,n);

%%

N = 100000000;
X = randi(n);
Xhist = zeros(N,1);
xhist = zeros(n,1);

xhist(X) = xhist(X)+1;
Xhist(1) = X;

for i=2:N
    
    if rand > alpha % then teleport
        X = randi(n);
    else
        if randi(2) == 1
            % Then we interpret our state as the last state,
            % and draw a random initial state from the distribution
            myhist = [X, Xhist(randi(i-1))];
        else
            % Then we interpret our state as the initial state, and
            % draw a random last state from the distribution
            myhist = [Xhist(randi(i-1)), X];
        end
        % now draw the next state
        Pnext = P(:,myhist(1),myhist(2));
        % Customize on the uniform case
        [nextX,~,probX] = find(Pnext);
        X = nextX(randi(numel(nextX)));
    end
    
    xhist(X) = xhist(X)+1;
    Xhist(i) = X;
end
xvec = xhist./N

%% Tensor PageRank
xtpr = tensorpr3(R,alpha).solve();

%% Hey, this seems like it works!

%% Let's try it for a general problem
clear;
load ../../tensors/ling1.mat
R = P;
n = size(R,1);
P = reshape(R,n,n,n);
Ps = cumsum(P); % useful for fast sampling

N = 1000000;
Xhist = zeros(N,1);
xhist = ones(n,1);

Xhist(1:n) = 1:n;
X = 1;

for i=n+1:N
    myhist = Xhist(randi(i-1));
    [xn,~,pn] = find(Ps(:,X,myhist));
    %Xn = xn(randsample(numel(xn),1,true,pn));
    Xn = xn(find(pn > rand,1,'first'));
    
    X = Xn;
    xhist(X) = xhist(X)+1;
    Xhist(i) = X;
end
xvec = xhist./N

    
%%
xtpr = tensorpr3(R,1).shifted

%% What happens if we have a non-convergent process
load ../../tensors/mtm4.mat
R = R4_19;
alpha = 0.99;
n = size(R,1);
P = reshape(R,n,n,n);

N = 1000000;
X = randi(n);
Xhist = zeros(N,1);
xhist = ones(n,1);

xhist(X) = xhist(X)+1;
Xhist(1:n) = 1:n;
Xhist(n+1) = X;

for i=n+2:N
    
    if rand > alpha % then teleport
        X = randi(n);
    else
        myhist = Xhist(randi(i-1));
        % now draw the next state
        Pnext = P(:,X,myhist);
        % Customize on the uniform case
        [nextX,~,probX] = find(Pnext);
        X = nextX(randi(numel(nextX)));
    end
    
    xhist(X) = xhist(X)+1;
    Xhist(i) = X;
end
xvec = xhist./N

%%
xtpr = tensorpr3(R,alpha).solve();
%%
exact = R4_Properties.R4_19.alpha99.sols
%%
% These were not the same!
% Check what we got from the power method
[~,~,~,xhist] = tensorpr3(R,alpha).power();
%%
mean(xhist(:,end-1000:end),2)
%%
% That also wasn't the same. What about the ODE?
x = ones(n,1)/n;
N = 100000;
xhist = zeros(n,N);
xhist(:,1) = x;
for i=2:N
    M = eye(n) - alpha*R*kron(eye(n),x);
    xpr = M\((1-alpha)*ones(n,1)/n);
    
    x = x +  (xpr - x);
    x = x/sum(x);
    xhist(:,i) = x;
end
%%
% Hmm... this seems to converge. What is going on? Let's try ODE45
fx = @(x) ((eye(n) - alpha*R*kron(eye(n),x))\((1-alpha)*ones(n,1)/n))-x;
[odet,odex] = ode45(@(t,x) fx(x) - sum(fx(x))*ones(n,1)/n, [0,100],ones(n,1)/n);
plot(odex)

%% 
% That shows we'll converge... 
% except how long might it take?

% Hmm... this seems to converge. What is going on? Let's try ODE45
fx = @(x) ((eye(n) - alpha*R*kron(eye(n),x))\((1-alpha)*ones(n,1)/n))-x;
[odet,odex] = ode45(@(t,x) 1/t*(fx(x) - sum(fx(x))*ones(n,1)/n), [1,1000000000],ones(n,1)/n,odeset('RelTol',1e-6,'AbsTol',1e-8));
plot(odet,odex);
set(gca,'XScale','log');

%%
% Let's see if a slightly forgetful stochastic process will do.


N = 1e10;
Nhist = 100000000;
X = randi(n);

Xhist = zeros(Nhist,1);
xhist = ones(n,1);

xhist(X) = xhist(X)+1;
Xhist(1:n) = 1:n;
Xhist(n+1) = X;

for i=n+2:N
    if mod(i,Nhist/10)==0
        fprintf('%g entries\n', i);
    end
    if rand > alpha % then teleport
        X = randi(n);
    else
        if i>Nhist
            myhist = Xhist(randi(Nhist));
        else
            Xhist(randi(i-1));
        end
        % now draw the next state
        Pnext = P(:,X,myhist);
        % Customize on the uniform case
        [nextX,~,probX] = find(Pnext);
        X = nextX(randi(numel(nextX)));
    end
    
    xhist(X) = xhist(X)+1;
    Xhist(mod(i,Nhist)+1) = X;
end
xvec = xhist./N