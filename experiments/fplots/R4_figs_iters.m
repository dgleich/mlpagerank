function R4_figs_iters
%% Load data
mats = load('../../tensors/mtm4.mat');
alpha = 0.99;
gamma = 1/2;
maxiter = 150;
n = 4;
v = ones(n,1)/n;
niter = 20000;
subset = mats.R4_mats;
subset = {'R4_12','R4_17','R4_18','R4_19'};
%%
for i=1:numel(subset)
    %%
    mat = subset{i};
    R = mats.(mat);
        
    tpr = tensorpr3(R,alpha,v);
    
    [~,~,flag1,xhist1] = tpr.shifted(0,'maxiter',20000);
    [~,~,flagg,xhistg] = tpr.shifted(gamma,'maxiter',20000);
    [~,~,flagi,xhisti] = tpr.inverseiter('maxiter',2000);
    [~,~,flagn,xhistn] = tpr.newton('maxiter',2000);
    [~,~,flagio,xhistio] = tpr.innout('maxiter',2000);
    
    clf;
    h = plot([1 1],[1 1],'-',[1 1],[1 1],'-',[1 1],[1 1],'-',[1 1],[1 1],'-',...
        'LineWidth',1.5);
    axis off;
    axis tight;
    legtext = cellstr(num2str((1:n)'));
    legend(legtext{:},'Location','WestOutside');
    legend boxoff;
    set(h,'Visible','off')
    axis tight;
    axis square;
    set_figure_size([2.5,0.5]);
    print(gcf,'R4_legend.eps','-depsc2');
    
    iterplot(mat,[v xhist1],maxiter,flag1,0);
    print(gcf,sprintf('%s-fixed.eps',mat),'-depsc2'); 
    iterplot(mat,[v xhistg],maxiter,flagg,0);
    print(gcf,sprintf('%s-shifted.eps',mat),'-depsc2');
    iterplot(mat,[v xhisti],maxiter,flagi,0);
    print(gcf,sprintf('%s-inverse.eps',mat),'-depsc2');
    iterplot(mat,[v xhistn],maxiter,flagn,0);
    print(gcf,sprintf('%s-newton.eps',mat),'-depsc2');
    iterplot(mat,[v xhistio],maxiter,flagio,0);
    print(gcf,sprintf('%s-innout.eps',mat),'-depsc2');
end