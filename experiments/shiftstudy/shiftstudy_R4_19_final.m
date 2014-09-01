%% Check if R4_19 converges for any shift.
% This will 

load ../../tensors/mtm4

R = R4_19;
n = size(R,1);
alpha = 0.99;
tpr = tensorpr3(R,alpha);

%% Final figure
clf; hold all; set(gca,'YScale','log'); set(gca,'XScale','log');
gammas = fliplr([0.50 0.553 0.554 0.5545 0.555 0.56 0.65 1]);
ng = numel(gammas);
niter = 1e6;
nskip = 200;
Results = {};
for i=1:ng
    [x,hist,flag] = tpr.shifted(gammas(i),'maxiter',niter);
    lasthist = [numel(hist) hist(end)];
    hist = hist(1:nskip*floor(numel(hist)/nskip));
    meanhist = mean(reshape(hist,nskip,numel(hist)/nskip));
    xmeanhist = nskip/2:nskip:numel(hist); % x coordinates of mean hist
    meanhist(end+1) = lasthist(2);
    xmeanhist(end+1) = lasthist(1);
    Results{i} = {[1:(nskip/2-1) xmeanhist], [hist(1:(nskip/2-1))' meanhist]};
    i
end
%%
save 'results_R4_19.mat' niter nskip gammas ng R alpha
%%
load results_R4_19
clf; hold all; set(gca,'YScale','log'); set(gca,'XScale','log');
for i=1:ng
    h=plot(Results{i}{1},Results{i}{2});
    ind = find(Results{i}{2} < 1e-5,1,'first');
    if ~isempty(ind)
        xlab = Results{i}{1}(ind);
        ylab = Results{i}{2}(ind);
    else
        xlab = Results{i}{1}(end);
        ylab = Results{i}{2}(end);
    end
    
    text(xlab*1.05,ylab*0.95,...
        sprintf('\\gamma=%g',gammas(i)),...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom',...
        'Color',get(h,'Color'));
end
    
pos = get(gca,'Position'); pos(3) = pos(3)*0.95; pos(2) = pos(2)+0.05;
set(gca,'Position',pos);
ylim([1e-6,1]);
xlabel('Residual');
box off;
ylabel('Iteration');


set_figure_size([5,3]);
print(gcf,'shiftstudy_R4_19.eps','-depsc2');
