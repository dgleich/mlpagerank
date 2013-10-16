function plotKappa(name,R,alpha,v)
I = eye(size(R,1));
step = @(x) inv(I - alpha/2.*R*(kron(I, x) + kron(x, I)))*v*(1-alpha);

%% Compute the solution trajectories.
[xtrue,xhist] = tensorRankNonShift(alpha,R,v,v);

%% Setup the simplex sampling.
[T,S,X,C] = simplex_samples(100);

%% Compute the solution trajectories on the simplex
xhist = [v xhist];
shist = C*xhist(1:2,:) + X(1:2,3)*ones(1,length(xhist));


%% Plot violation of the norm bound on the simple
% compute the violation of our norm bound
fxp = @(x,s) norm(step(s) - s,1)/norm(s-x,1);
fx = @(x) fxp(x,step(x));
nj = [];
for i=1:size(S,2);
    nj(i) = fx(S(:,i));
end
clf; hold all;
scolor = get(gcf,'Color');
scolor = (scolor + [1,1,1])/2;
%patch(X(1,:),X(2,:),scolor,'EdgeColor',scolor);
scatter(T(1,:),T(2,:),25,nj,'filled');
cmapsetup;
colorbar;

hold on;
plot(shist(1,:),shist(2,:),'g-');
hold off;


axis off;
axis tight;
axis square;
axis equal;

title(sprintf('max = %f\n', max(nj)), 'FontSize', 15);
set(gca, 'FontSize', 15);
set_figure_size([3 3]);
print(gcf,sprintf('%s-stepwise-bound.eps',name), '-depsc2');



end