function tensorposter_study_figures(name,R,alpha,v,gamma)
step = @(x) alpha*R*kron(x,x) + (1-alpha)*v;

%% Compute the solution trajectories.
[xtrue1,xhist1] = tensorRank(alpha,R,v,v,1);
[xtrueg,xhistg] = tensorRank(alpha,R,v,v,gamma);

%% Setup the simplex sampling.
[T,S,X,C] = simplex_samples(100);

%% Compute the solution trajectories on the simplex
xhistg = [v xhistg];
shistg = C*xhistg(1:2,:) + X(1:2,3)*ones(1,length(xhistg));

xhist1 = [v xhist1];
shist1 = C*xhist1(1:2,:) + X(1:2,3)*ones(1,length(xhist1));


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
plot(shist1(1,:),shist1(2,:),'g-');
plot(shistg(1,:),shistg(2,:),'Color',[1,0.4,0.6]);
hold off;


axis off;
axis tight;
axis square;
axis equal;

title(sprintf('max = %f\n', max(nj)));

set_figure_size([4 4]);
print(gcf,sprintf('%s-stepwise-bound.eps',name), '-depsc2');

%% Plot the second largest eigenvalue magnitude on the simplex
fxp = @(s) s(2);
prodejac = @(x, alpha, R, v) alpha*R*(kron(x, eye(size(x,1))) + kron(eye(size(x,1)), x)) - eye(size(x,1)); 
%fx = @(x) fxp(sort(abs(eig(prodejac(x,alpha,R,v))),'descend'));
% second largest eigenvalue, not the magnitude
fx = @(x) fxp(sort(eig(prodejac(x,alpha,R,v)),'descend'));
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

hold on;
plot(shist1(1,:),shist1(2,:),'g-');
plot(shistg(1,:),shistg(2,:),'Color',[1,0.4,0.6]);
hold off;

axis off;
axis tight;
axis square;
axis equal;

title(sprintf('max = %f\n', max(nj)));

set_figure_size([4 4]);
print(gcf,sprintf('%s-jac.eps',name), '-depsc2');

