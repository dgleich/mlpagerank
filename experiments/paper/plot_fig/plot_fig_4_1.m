%%
addpath('../../../');
alpha = 0.99;
v = ones(3,1)./3;
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
[~, xhist] = tensorRank(alpha, R1, v, v, 0);
xhist = [v xhist]';
plot(xhist(1:50, :));
set(gca, 'FontSize', 15);
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3');
axis([1 50 0 0.7]);
set_figure_size([3 3]);
name = 'Richardson';
print(gcf,sprintf('%s.eps',name), '-depsc2');
