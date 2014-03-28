function smap = tetrahedron_plot_kappa(name,tpr)
%% Setup the simplex sampling.
[T,S,X,C] = tetrahedron_samples(50);
% S are points on the simplex
% X are points projected to the tetrahedron
% C helps project
%% 
smap = @(x) C*x(1:3, :) + T(:, 4)*ones(1, size(x, 2));
%% compute the function value for the volume
step = @(x) tpr.residual(x) + x;
% compute the violation of our norm bound
fxp = @(x,s) norm(step(s) - s,1)/norm(s-x,1);
fx = @(x) fxp(x,step(x));
nj = zeros(size(S,2),1); % at some point nj was the norm of the Jacobian
for i = 1:size(S,2)
    nj(i) = fx(S(:,i));
end
%% kappa value
clf;
hold all;

ind = find(nj>1);
P = X(:,ind); % points to plot
vnj = nj(ind); % values nj
hold on;
plot3([0 0 1 0 1 1 1 0], [0 1 1 0 0 1 0 1], [1 0 1 1 0 1 0 0], 'k--');
scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
cmapsetup;
colorbar;

axis off;
axis tight;
axis square;
axis equal;

title(sprintf('%s max = %f\n', name, max(nj)), 'FontSize', 12, ...
    'Interpreter','none');
set(gca, 'FontSize', 12);


set_figure_size([3 3]);

