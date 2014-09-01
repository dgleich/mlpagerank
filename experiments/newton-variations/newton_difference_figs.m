%% Make a figure of the difference between projected newton and otherwise
M = load('../../tensors/mtm4');
alpha = 0.9;
R = M.R4_11;
tpr = tensorpr3(R,alpha);
%%
[~,~,~,xhist] = newton_nonorm(tpr,'x0',0);
clf;
plot(xhist','.-','MarkerSize',8)
ylim([0,0.1]);
xlim([1,5]);
xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;
print(gcf,'vanilla-newton.eps','-depsc2');
%%
[~,~,~,xhist] = newton_project(tpr,'x0',0);
clf;
plot(xhist','.-','MarkerSize',8)
ylim([0,1]);
xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;
print(gcf,'stochastic-newton.eps','-depsc2');

%%
!cp *.eps ~/Dropbox/publications/tensorpr-shared/figures/