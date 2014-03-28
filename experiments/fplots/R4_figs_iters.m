function R4_figs_iters
%% Load data
mats = load('../../tensors/mtm4.mat');
alpha = 0.99;
gamma = 2/3;
maxiter = 150;
n = 4;
v = ones(n,1)/n;
%%
for i=1:numel(mats.R4_mats)
    %%
    mat = mats.R4_mats{i};
    R = mats.(mat);
        
    tpr = tensorpr3(R,alpha,v);
    
    [~,~,flag1,xhist1] = tpr.shifted(0);
    [~,~,flagg,xhistg] = tpr.shifted(gamma);
    [~,~,flagi,xhisti] = tpr.inverseiter;
    [~,~,flagn,xhistn] = tpr.newton;
    [~,~,flagio,xhistio] = tpr.innout;
    
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