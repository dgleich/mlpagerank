function R3_figs
%% Load data
mats = load('../../tensors/mtm3.mat');
alpha = 0.99;
gamma = 1/2;
v = ones(3,1)/3;
maxiter = 100;
n = 3;
%%
for i=1:numel(mats.R3_mats)
    %%
    mat = mats.R3_mats{i};
    R = mats.(mat);
    
    [~,~,flag1,xhist1] = tensorpr3(R,alpha,v).shifted(0);
    [~,~,flagg,xhistg] = tensorpr3(R,alpha,v).shifted(gamma);
    [~,~,flagi,xhisti] = tensorpr3(R,alpha,v).inverseiter;
    [~,~,flagn,xhistn] = tensorpr3(R,alpha,v).newton;
    [~,~,flagio,xhistio] = tensorpr3(R,alpha,v).innout;
    
    clf;
    h = plot([1 1],[1 1],'-',[1 1],[1 1],'-',[1 1],[1 1],'-','LineWidth',1.5);
    axis off;
    axis tight;
    legtext = cellstr(num2str((1:n)'));
    legend(legtext{:},'Location','WestOutside');
    legend boxoff;
    set(h,'Visible','off')
    axis tight;
    axis square;
    set_figure_size([2.5,0.5]);
    print(gcf,'R3_legend.eps','-depsc2');
    
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

    smap = triangle_plot_lambda2(mat,R,alpha,v);
    myplot = @(x,varargin) plot(x(1,:),x(2,:),varargin{:},'LineWidth',1.5);
    hold all;
    myplot(smap([v xhist1]),'g-');
    myplot(smap([v xhistg]),'r-','Color',[1,0.4,0.6]);
    hold off;
    print(gcf,sprintf('%s-jac.eps',mat),'-depsc2');
    
    smap = triangle_plot_kappa(mat,R,alpha,v);
    hold all;
    myplot(smap([v xhist1]),'g-');
    myplot(smap([v xhistg]),'r-','Color',[1,0.4,0.6]);
    hold off;
    print(gcf,sprintf('%s-kappa.eps',mat),'-depsc2');
    
    %triangle_plot_kappa(mat,R,0.99,ones(3,1)/3,1/2);
end