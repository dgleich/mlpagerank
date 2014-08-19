%% Check if there is any corresponds with an edge reinforced type process
% and our crazy stochastic process

% We are going to run an edge reinforced random walk on the graph of the
% second-order chain to see what happens with it.
P = zeros(2,2,2);
P(1,1,1) = 1;
P(2,1,1) = 1;
P(2,2,1) = 1;
P(1,2,1) = 1;
P(1,2,2) = 1;
P(1,1,2) = 1;

R = reshape(P,2,4)

M = tensor2markov(R)


%%
H = M;
X = 1;
N = 100000;
States = ones(4,1);
for i=1:N
    States(X) = States(X)+1;
    [xn,~,pn] = find(H(:,X));
    Xn = xn(randsample(numel(xn),1,true,pn));
    H(Xn,X) = H(Xn,X) + 1;
    X = Xn;
end
H    
    
%% Now let's run our crazy process and track what it does in the tensor 
% product space
Hc = P;

N = 100000;
X = 1; % we are in state 1


Xhist = zeros(N+2,1);
xhist = ones(2,1);
Xhist(1:2) = [2,1];


for i=3:N
    myhist = Xhist(randi(i-1));
    Pnext = P(:,X,myhist);
    [xn,~,pn] = find(P(:,X,myhist));
    Xn = xn(randsample(numel(xn),1,true,pn));
    % we are taking a transition from X->Xn, last was Xhist(i-1)
    %Hc(Xn,X,Xhist(i-2)) = Hc(Xn,X,Xhist(i-2)) + 1;
    Hc(Xn,X,myhist) = Hc(Xn,X,myhist) + 1;
    
    X = Xn;
    
    xhist(X) = xhist(X)+1;
    Xhist(i) = X;
end
xvec = xhist./N

Hcm = tensor2markov(reshape(Hc,2,4))