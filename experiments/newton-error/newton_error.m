%% Test a conjecture about the error in Newton's method
load tensors/example1
alpha = 0.4;
n = size(R,1);
v = zeros(n,1);
S = tensorpr3(R,alpha);
x = S.solve('tol',1e-16);

%% Compute one step of Newton's method
x0 = ones(n,1)./n;
x1 = x0 - (S.jacobian(x0)-eye(n))\(S.residual(x0));

%% Check our prediction of the error
% Based on our analysis, we should have:
%   error at k+1 = (S.jacobian(x at k)-eye(n)) \ 
%                          (alpha*R*[(error at k) kron (error at k)]
% but this could be wrong. 
e0 = x - x0; % error at first ste
e1 = x - x1; %

ep = (S.jacobian(x0) - eye(n)) \ (alpha*R*kron(e0,e0))
%%
% Hmm... it seems we are off by a negative sign
ep = -(S.jacobian(x0) - eye(n)) \ (alpha*R*kron(e0,e0))
e1 - ep

%% Let's check at each iteration now.
xcur = ones(n,1)./n;
e = x - xcur;
for i=1:10
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    norm(e-ep)
end

%% That holds for all iterations. Awesome. 
xcur = ones(n,1)./n;
e = x - xcur;
for i=1:10
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    [norm(inv((S.jacobian(xiter) - eye(n))),1) 
        norm(alpha*R*kron(eold,eold),1) ]
end    

%% Now check Newton's method with a maximum initial error to see what happens

load tensors/mtm3.mat
R = R3_4;
alpha = 0.49;
n = size(R,1);
v = zeros(n,1);
v(3) = 1;
S = tensorpr3(R,alpha,v);
x = S.solve('tol',1e-16);

%% That holds for all iterations. Awesome. 
xcur = zeros(n,1);
xcur(1) = 1;
e = x - xcur;
for i=1:1
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    [norm(e); 
        norm(alpha*R*kron(eold,eold),1); norm(e-ep)]
end    

%% Test if we keep things positive if x0 = 0.

load tensors/mtm3.mat
R = R3_4;
alpha = 0.49;
n = size(R,1);
v = zeros(n,1);
v(3) = 1;
S = tensorpr3(R,alpha,v);
x = S.solve('tol',1e-16);

xcur = zeros(n,1);
e = x - xcur;
fprintf('\n\n\n\nStart\n\n\n\n');
for i=1:10
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    [sum(e); norm(e); sum(xcur - xiter); 
        norm(alpha*R*kron(eold,eold),1); norm(xcur - xiter,1)]'
end    

% Awesome, we do, because sum(xcur-xiter) = norm(xcur-xiter)

%% Look at the sum of the changes in pk = x_k+1 - x_k
load tensors/mtm3.mat
R = R3_4;
alpha = 0.49;
n = size(R,1);
v = zeros(n,1);
v(3) = 1;
S = tensorpr3(R,alpha,v);
x = S.solve('tol',1e-16);

xcur = zeros(n,1);
e = x - xcur;
fprintf('\n\n\n\nStart\n\n\n\n');
for i=1:10
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    [sum(S.residual(xcur)) alpha*sum(xcur-xiter)^2 sum(xcur-xiter) sum(S.residual(xiter))/(1-2*sum(xiter)*alpha) ]'
end    

%% Now, we develop recurrences for each of the terms
load tensors/mtm3.mat
R = R3_4;
alpha = 0.49;
n = size(R,1);
v = zeros(n,1);
v(3) = 1;
S = tensorpr3(R,alpha,v);
x = S.solve('tol',1e-16);

xcur = zeros(n,1);
e = x - xcur;
fprintf('\n\n\n\nStart\n\n\n\n');
for i=1:3
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    fk = sum(S.residual(xiter));
    fk1 = sum(S.residual(xcur));
    sumx = (1-sqrt(1-4*alpha*(1-alpha - fk)))/(2*alpha);
    fk12 = alpha*fk^2/((1-2*alpha)^2 + 4*alpha*fk);
    [sum(xcur) sumx alpha*(fk/(1-2*alpha*sumx))^2 fk1 fk12 ]
end    
%% The recurrence
%a = (1/2 - 0.25);
a = 1/2 - 0.4;
f = (1-a);
fk = [];
for i=1:25
    f = a*f^2/((1-2*a)^2 + 4*a*f);
    fk(i) = f;
    s = '';
    %gk = (1-a)^2*(2*a)^(2^i);
    theta = a/((1-2*a)^2);
    %c = (a*theta - 1)/(2*theta);
    c = max(1-2*a,2*a);
    %c = max(2*a);
    gk = min(theta*(c)^(2^i),(1/4)^(i-1)*(1-a)^2*a);
    if gk < f, s = '****'; end
    fprintf('%4i  %.20f  %.20f  %.20f %s\n', i, f, gk, (1/4)^(i-1)*(1-a)^2*a, s);
    %fprintf('%4i  %.20f  %.20f\n', i, f, (1-2*a)^2/(a)*((1-2*a)^(2^(i))));
end

% Try and determine the rate of convergence
rates = fk(2:end) ./ fk(1:end-1).^2;
rates = rates(isfinite(rates));
fprintf('rate = %.6f   a = %.6f   guess = %.6f\n', rates(end), a, theta);


%% The simple convergence result
% Yongyang pointed out that there is a simple result that: f_{k+1} <=
% 1/4*f_k, which gives linear convergence. There is probably a quadratic
% result too, but linear will do for now.

%% Let's justify the result of using the negative root in the proof.
% At one step of the recurrence, I experimentally confirmed that we use the
% negative root of a quadratic, but now we need to prove that.

load tensors/mtm3.mat
R = R3_4;
alpha = 0.49;
n = size(R,1);
v = zeros(n,1);
v(3) = 1;
S = tensorpr3(R,alpha,v);

xcur = zeros(n,1);
z = 0;
e = x - xcur;
fprintf('\n\n\n\nStart\n\n\n\n');
for i=1:10
    xiter = xcur;
    eold = e;
    xcur = xcur - (S.jacobian(xcur)-eye(n))\(S.residual(xcur));
    e = x - xcur;
    ep = -(S.jacobian(xiter) - eye(n)) \ (alpha*R*kron(eold,eold));
    fk = sum(S.residual(xiter));
    fk1 = sum(S.residual(xcur));
    sumx = (1-sqrt(1-4*alpha*(1-alpha - fk)))/(2*alpha);
    sumx2 = (1+sqrt(1-4*alpha*(1-alpha - fk)))/(2*alpha);
    z = z + 1/(1-2*alpha*z)*fk;
    [sum(xcur) sumx sumx2 z ]
end 

%% Okay, test our final convergence bound
for ai=linspace(0,0.5,10000)
    a = ai;
    f = (1-a);
    fk = [];
    for i=1:25
        f = a*f^2/((1-2*a)^2 + 4*a*f);
        fk(i) = f;
        s = '';
        %gk = (1-a)^2*(2*a)^(2^i);
        theta = a/((1-2*a)^2);
        %c = (a*theta - 1)/(2*theta);
        c = max(1-2*a,2*a);
        gk = min(theta*(c)^(2^i),(1/4)^(i-1)*(1-a)^2*a);
        if gk < (1-sqrt(eps))*f, 
            s = '****'; 
            fprintf('%.3f  %4i  %.20f  %.20f  %.20f %s\n', a, i, f, gk, (1/4)^(i-1)*(1-a)^2*a, s);
        end
        %fprintf('%4i  %.20f  %.20f\n', i, f, (1-2*a)^2/(a)*((1-2*a)^(2^(i))));
    end
end

%% Try a relationship based on the quadratic equation
a = 0.49;
f = (1-a);
fk = [];
for i=1:25
    finit = f;
    f = a*f^2/((1-2*a)^2 + 4*a*f);
    fk(i) = f;
    s = '';
    %gk = (1-a)^2*(2*a)^(2^i);
    theta = a/((1-2*a)^2);
    %c = (a*theta - 1)/(2*theta);
    c = max(1-2*a,2*a);
    %gk = min(theta*(c)^(2^i),(1/4)^(i-1)*(1-a)^2*a);
    gk = (-(1-2*a)^2 - sqrt((1-2*a)^4 + 16*a^2*finit^2))/(8*a);
    if gk < (1-sqrt(eps))*f, 
        s = '****'; 
    else
        s = '';
    end
    fprintf('%.3f  %4i  %.20f  %.20f  %.20f %s\n', a, i, f, gk, 1, s);
   
    %fprintf('%4i  %.20f  %.20f\n', i, f, (1-2*a)^2/(a)*((1-2*a)^(2^(i))));
end
