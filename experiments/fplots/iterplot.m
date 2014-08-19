function h = iterplot(name,iters,maxiter,flag,printlegend)

xmax = maxiter;
if maxiter > size(iters,2)
    maxiter = size(iters,2);
end
n = size(iters,1);
plot(iters(:,1:maxiter)','LineWidth',1.25);
xlim([1,xmax]);
set(gca,'LineWidth',0.85);
set(gca,'FontSize',10);
ylim([0,1]);
ctext = 'Converged';
if ~flag
    ctext = 'Did not converge';
end
if maxiter < xmax
    text(xmax,1,ctext,'HorizontalAlignment','right','VerticalAlignment','top');
else
    text(1,1,ctext,'HorizontalAlignment','left','VerticalAlignment','top');
end
if printlegend
    legtext = cellstr(num2str((1:n)'));
    legend(legtext{:},'Location','OutsideWest');
end
box off;
set_figure_size([2.5,2.5]);
