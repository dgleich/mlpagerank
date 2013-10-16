function plotJacobianEig(name, R, alpha, v)
% This function is used to generate the triangle plot
% for the tensorRankNonShift method.
% Specifically, the plot outputs the second largest eigenvalue
% values on the simplex.	
I = eye(size(v, 1));

%% Compute the soln trajactories
[xtrue, xhist] = tensorRankNonShift(alpha, R, v, v);

%% Setup the simplex sampling
[T, S, X, C] = simplex_samples(100);

%% Compute the soln trajectories on the simplex
xhist = [v xhist];
shist = C*xhist(1:2, :) + X(1:2, 3)*ones(1, length(xhist));

%% plot the second largest eigenvalue on the simplex
fxp = @(s) s(2);
prodejac = @(x, alpha, R, v) getJ(alpha, R, v, x);
fx = @(x) fxp(sort(eig(prodejac(x, alpha, R, v)), 'descend'));
nj = [];
for i=1:size(S,2)
    nj(i) = fx(S(:,i));
end
clf; hold all;
scolor = get(gcf, 'Color');
scolor = (scolor+[1,1,1])/2;
patch(X(1,:), X(2,:),scolor, 'EdgeColor', scolor);
scatter(T(1,:), T(2, :), 25, nj, 'filled');
cmapsetup;
caxis([-1.5 0.75]);
hcb = colorbar;
set(hcb, 'YTick', [-1.5 -0.75 0 0.75]);

hold on;
plot(shist(1,:), shist(2,:), 'g-');
hold off;

axis off;
axis tight;
axis square;
axis equal;

title(sprintf('max = %f\n', max(nj)), 'FontSize', 15);

set_figure_size([3 3]);
print(gcf, sprintf('%s-jac.eps', name), '-depsc2');
end