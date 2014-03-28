%% Test tetrahedron samples
[T,S,X,C] = tetrahedron_samples(50);
% S are points on the simplex
% X are points projected to the tetrahedron
% C helps project
%% 
smap = @(x) C*x(1:3, :) + T(:, 4)*ones(1, size(x, 2));

%% compute the function value for the volume
P = X; % points to plot
vnj = max(S);
thresh = 0.8;
P = X(:,vnj >thresh);
vnj = vnj(vnj > thresh);
clf; hold on;
scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
h = trimesh([1 2 3; 1 2 4; 1 3 4; 2 3 4], T(1,:), T(2,:), T(3,:),[0,0,0,0]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
%scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
cmapsetup_large

%% compute the function value for the volume
P = X; % points to plot
vnj = sum(smap(S).*smap(S));
thresh = 0.4;
P = X(:,vnj >thresh);
vnj = vnj(vnj > thresh);
clf; hold on;
scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
h = trimesh([1 2 3; 1 2 4; 1 3 4; 2 3 4], T(1,:), T(2,:), T(3,:),[0,0,0,0]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
axis off;
axis tight;
axis equal;
view([29,68]);
camzoom(1.5);
%%
camlight headlight;
lighting gouraud

%scatter3(P(1,:), P(2,:), P(3,:), 25, vnj, 'filled');
cmapsetup_large

