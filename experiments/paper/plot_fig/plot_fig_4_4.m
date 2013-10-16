addpath('../../../');
alpha = 0.99;
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
%% fig 4.4(a)
v = ones(3,1)./3;
[~, xhist] = tensorRankNonShift(alpha, R1, v, v);
xhist = [v xhist]';
plot(xhist);
set(gca, 'FontSize', 15);
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3');
axis([1 28 0 0.58]);
set_figure_size([3 3]);
name = 'Non-shift-converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');
%% fig 4.4(b)
addpath('../../../test_case');
load('6x6x6_not_converge_non_shift.mat', 'R1');
R2 = R1;
v = ones(6,1)./6;
[~, xhist] = tensorRankNonShift(alpha, R2, v, v);
xhist = [v xhist]';
plot(xhist(1:100, :));
set(gca, 'FontSize', 15);
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6');
axis([1 100 0 0.71]);
set_figure_size([3 3]);
name = 'Non-shift-non-converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');
%% fig 4.4(c)
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
name = 'Non-shift';
v = ones(3,1)./3;
plotKappa(name, R1, alpha, v);
%% fig 4.4(d)
name = 'Non-shift';
v = ones(3,1) ./ 3;
plotJacobianEig(name, R1, alpha, v);