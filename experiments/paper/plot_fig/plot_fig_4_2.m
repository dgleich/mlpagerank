%% fig 4.2(a)
addpath('../../../');
alpha = 0.99;
v = ones(3,1)./3;
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
[~, xhist] = tensorRank(alpha, R1, v, v);
xhist = [v xhist]';
plot(xhist(1:50, :));
set(gca, 'FontSize', 15);
title('shift \gamma = 1/2');
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3');
axis([1 50 0 0.58]);
set_figure_size([3 3]);
name = 'shift_converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');
%% fig 4.2(b)
alpha = 0.99;
v = ones(6,1)./6;
addpath('../../../test_case');
load('6x6x6_not_converge_non_shift.mat', 'R1');
R2 = R1;
[~, xhist] = tensorRank(alpha, R2, v, v);
xhist = [v xhist]';
plot(xhist(1:200, :));
set(gca, 'FontSize', 15);
title('shift \gamma = 1/2');
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6');
axis([1 200 0 0.8]);
set_figure_size([3 3]);
name = 'shift_non_converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');