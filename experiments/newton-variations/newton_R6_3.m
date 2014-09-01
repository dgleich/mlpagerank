%% Newton's method shows some really interesting behavior on R6_3.
R6 = load('../../tensors/mtm6.mat');
R = R6.R6_3;
alpha = 0.99;
tpr = tensorpr3(R,0.99);
gamma = R6.R6_Properties.R6_3.alpha99.gamma;

%% True solution
Y = [   0.043820721946272
   0.002224192630620
   0.009256490884022
   0.819168263512464
   0.031217440669761
   0.094312890356862];


%%
[~,hist] = newton_project(tpr);
loglog(hist)
%%
[~,hist] = newton_prox(tpr);
loglog(hist)
%%
[~,hist] = newton_project(tpr);
clf;
semilogy(hist(1:20),'.-'); hold all;
plot(8,hist(8),'ro','MarkerSize',6);
ylabel('Residual');
xlabel('Iteration');
set_figure_size([2,2]);
ylim([1e-8,1]);
set(gca,'YTick',[1e-8,1e-6,1e-4,1e-2,1]);
box off;
print(gcf,'newton_failure_R6_3.eps','-depsc2');
%%
[~,hist,~,xhist] = newton_project(tpr,'maxiter',8);
hist(end)
x0 = xhist(:,end);
norm(tpr.residual(x0),1)
eig(tpr.jacobian(x0)-eye(6))

%% Try other algorithms starting at x0
[~,hist,~,xhist] = tpr.innout('x0',x0);
loglog(hist)

%% 


%% Solve for the step without that eigenvalue
J = tpr.jacobian(x0) - eye(6);
[V,D] = eig(J);
D(2,2) = 0;
J2 = V*D*inv(V)
xn = J2 \ tpr.residual(x0);

%% We can find the true solution by random sampling
for N=1:100
    [x,hist,flag,xhist] = tpr.newton('randinit',true);
    if flag
        fprintf('Found!\n');
        xhist(:,1)
        break;
    end
end
%%
Xset = xhist(:,end-10:end);
for i=1:size(Xset,2)
    [tpr.residual(Xset(:,i)); eig(tpr.jacobian(Xset(:,i))-eye(6))]
end