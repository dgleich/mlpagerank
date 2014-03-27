function smap = triangle_plot_kappa(name,R,alpha,v)
%% Setup the simplex sampling.
[T,S,X,C] = simplex_samples(100);

%% Compute the solution trajectories on the simplex
smap = @(xhist) C*xhist(1:2,:) + X(1:2,3)*ones(1,length(xhist));


%% Plot violation of the norm bound on the simple
step = @(x) alpha*R*kron(x,x) + (1-alpha)*v;
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
patch(X(1,:),X(2,:),scolor,'EdgeColor',scolor);
scatter(T(1,:),T(2,:),25,nj,'filled');
cmapsetup;
colorbar;

%hold on;
%plot(shist1(1,:),shist1(2,:),'g-');
%plot(shistg(1,:),shistg(2,:),'Color',[1,0.4,0.6]);
%hold off;

axis off;
axis tight;
axis square;
axis equal;

title(sprintf('%s max = %f\n', name, max(nj)), 'FontSize', 12, ...
    'Interpreter','none');
set(gca, 'FontSize', 12);
set_figure_size([3 3]);
%print(gcf,sprintf('%s-stepwise-bound.eps',name), '-depsc2');