function h = iterplot(name,iters,maxiter,printlegend)

xmax = maxiter;
if maxiter > size(iters,2)
    maxiter = size(iters,2);
end
n = size(iters,1);
plot(iters(:,1:maxiter)','LineWidth',1.5);
xlim([1,xmax]);
set(gca,'LineWidth',0.85);
set(gca,'FontSize',10);
ylim([0,1]);
if printlegend
    legtext = cellstr(num2str((1:n)'));
    legend(legtext{:},'Location','OutsideWest');
end
box off;
set_figure_size([2.5,2.5]);
