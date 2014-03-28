function smap = tetrahedron_plot_lambda2(name,tpr)
%% Setup the simplex sampling.
[T,S,X,C] = tetrahedron_samples(50);
% S are points on the simplex
% X are points projected to the tetrahedron
% C helps project
%% 
smap = @(x) C*x(1:3, :) + T(:, 4)*ones(1, size(x, 2));
%% compute the function value for the volume
fxp = @(s) s(2);
fx = @(x) fxp(sort(real(eig(tpr.jacobian(x)-eye(4))),'descend'));
nj = zeros(size(S,2),1); % at some point nj was the norm of the Jacobian
for i = 1:size(S,2)
    nj(i) = fx(S(:,i));
end
%% kappa value
clf;
hold all;

ind = find(nj>0);
P = X(:,ind); % points to plot
vnj = nj(ind); % values nj
hold on;
h = trimesh([1 2 3; 1 2 4; 1 3 4; 2 3 4], T(1,:), T(2,:), T(3,:),[0,0,0,0]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
hs = scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
cmapsetup_large;
%caxis([0 0.75]);
colorbar;

axis off;
axis tight;
axis equal;

view([-19.3592,-30.3143]);
camtarget([0,0,0]);
camzoom(2);


title(sprintf('%s max = %f\n', name, max(nj)), 'FontSize', 12, ...
    'Interpreter','none');
set(gca, 'FontSize', 12);


set_figure_size([3 3]);

