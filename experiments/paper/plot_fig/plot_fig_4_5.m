%% fig 4.5(a)
addpath('../../../');
alpha = 0.99;
v = ones(3,1)./3;
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
[~, xhist] = tensorNewton(alpha, R1, v, v);
xhist = [v xhist]';
plot(xhist);
set(gca, 'FontSize', 15);
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3');
axis([1 6 0 0.58]);
set_figure_size([3 3]);
name = 'Newton-converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');
%% fig 4.5(b)
addpath('../../../test_case/');
load('6x6x6_not_converge_non_shift.mat', 'R1');
R2 = R1;
v = ones(6,1)./6;
[~, xhist] = tensorNewton(alpha, R2, v, randv(6));
xhist = [v xhist]';
if max(size(xhist) > 200)
    plot(xhist(1:200, :));
else
    plot(xhist);
end
set(gca, 'FontSize', 15);
xlabel('Iteration(k)');
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6');
if max(size(xhist) > 200)
    axis([1 200 -8 6]);
end
set_figure_size([3 3]);
name = 'Newton-non-converge';
print(gcf,sprintf('%s.eps',name), '-depsc2');