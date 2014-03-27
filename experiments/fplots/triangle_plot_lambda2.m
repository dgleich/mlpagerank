function smap = triangle_plot_lambda2(name,R,alpha,v,varargin)
%% Setup the simplex sampling.
[T,S,X,C] = simplex_samples(100);

%% Compute the solution trajectories on the simplex
smap = @(xhist) C*xhist(1:2,:) + X(1:2,3)*ones(1,length(xhist));

%% Plot the second largest eigenvalue magnitude on the simplexfxp = @(s) s(2);
prodejac = @(x, alpha, R, v) alpha*R*(kron(x, eye(size(x,1))) + kron(eye(size(x,1)), x)) - eye(size(x,1)); 
%fx = @(x) fxp(sort(abs(eig(prodejac(x,alpha,R,v))),'descend'));
% second largest eigenvalue, not the magnitude
fxp = @(s) s(2);
fx = @(x) fxp(sort(real(eig(prodejac(x,alpha,R,v))),'descend'));
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
caxis([-1.5 0.75]);
hcb = colorbar;
set(hcb, 'YTick', [-1.5 -0.75 0 0.75]);

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
%print(gcf,sprintf('%s-jac.eps',name), '-depsc2');
end