%% This experiment makes figures to illustrate each algorithm.

load ../tensors/mtm3.mat
load ../tensors/mtm4.mat

%%
R = R3_1;
alpha = 0.95;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.power('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)

ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:3)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;
print('power-R3_1-95.eps','-depsc2');

%%
R = R3_1;
alpha = 0.96;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.power('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)

ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:3)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('power-R3_1-96.eps','-depsc2');

%% This same problem works with the shifted iteration
R = R3_1;
alpha = 0.96;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.shifted;

clf;
semilogx(xhist','.-','MarkerSize',8)

ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:3)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('shifted-R3_1-96.eps','-depsc2');

%% Actually all of the 3x3x3 problems converge with the shifted method

load ../tensors/mtm4

%% This problem does not converge with the shifted method
R = R4_11;
alpha = 0.97;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.shifted('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('shifted-R4_11-97.eps','-depsc2');


%% But it does with the innout
R = R4_11;
alpha = 0.97;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.innout('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('innout-R4_11-97.eps','-depsc2');

%% But then it doesn't again with the innout method
R = R4_11;
alpha = 0.99;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.innout('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('innout-R4_11-99.eps','-depsc2');

%% Inverse iteration on this same problem
R = R4_11;
alpha = 0.97;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.inverseiter('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('inverse-R4_11-97.eps','-depsc2');

%%
R = R4_11;
alpha = 0.99;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.inverseiter('maxiter',1000);

clf;
semilogx(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])
set(gca,'XTick',logticks(size(xhist,2)));

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('inverse-R4_11-99.eps','-depsc2');

%% Newton
R = R4_11;
alpha = 0.97;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.newton('maxiter',1000);

clf;
plot(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('newton-R4_11-97.eps','-depsc2');

%% Newton
R = R4_11;
alpha = 0.99;
tpr = tensorpr3(R,alpha);
[~,~,~,xhist] = tpr.newton('maxiter',1000);

clf;
plot(xhist','.-','MarkerSize',8)


ylim([0,1.5]);
set(gca,'YTick',[0:0.1:1]);
xlim([1,size(xhist,2)])

xlabel('Iteration');
ylabel('Solution');
legend(cellstr(num2str((1:4)')),'Location','Northwest');
legend boxoff;
set_figure_size([3,3]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);
for i=1:size(xhist,2)
    reshist(i) = norm(tpr.residual(xhist(:,i)),1);
end

semilogy(reshist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlabel('Iteration');
ylabel('Residual');
box off;

print('newton-R4_11-99.eps','-depsc2');

%% Print out tables
print_tensor3(spones(R3_1))
print_tensor3(spones(R4_11))