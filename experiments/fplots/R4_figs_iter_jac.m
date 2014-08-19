function R4_figs
%% Load data
mats = load('../../tensors/mtm4.mat');
alpha = 0.99;
gamma = 1/2;
maxiter = 100;
n = 4;
v = ones(n,1)/n;
subset = {'R4_12','R4_17','R4_18','R4_19'};
%%
for i=1:numel(subset)
    %%
    %mat = mats.R4_mats{i};
    mat = subset{i};
    R = mats.(mat);
        
    tpr = tensorpr3(R,alpha,v);
    %smap = tetrahedron_plot_kappa(mat,tpr);
    %tetrahedron_plot_lambda2(mat,tpr);
        
    %[~,~,~,xhist1] = tpr.shifted(0);
    [~,~,~,xhistg] = tpr.shifted(gamma);
    [~,~,~,xhisti] = tpr.inverseiter;
    [~,~,~,xhistn] = tpr.newton;
    [~,~,~,xhistio] = tpr.innout;
    
    fxp = @(s) s(2);
    fx = @(x) fxp(sort(real(eig(tpr.jacobian(x)-eye(4))),'descend'));
    evaljac = @(X) cellfun(fx,num2cell(X,1));
    
    clf; hold all;
    plot(evaljac(xhistg));
    plot(evaljac(xhisti));
    plot(evaljac(xhistn));
    plot(evaljac(xhistio));
    
    pause;
end