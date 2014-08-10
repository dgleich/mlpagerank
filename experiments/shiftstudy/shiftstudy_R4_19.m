%% Check if R4_19 converges for any shift.
% This will 

load ../../tensors/mtm4

R = R4_19;
n = size(R,1);
alpha = 0.99;
tpr = tensorpr3(R,alpha);
ng = 100;
gammas = linspace(0.5,1,ng+1);
gammas = gammas(2:end-1); % remove end points
convflag = zeros(ng-1,1);

%%
for i=1:ng-1
    [x,hist,flag] = tpr.shifted(gammas(i),'maxiter',20000);
    convflag(i) = flag;
    semilogy(hist); hold on;
    semilogy(100:200:numel(hist), mean(reshape(hist,200,numel(hist)/200)),...
        'r-','LineWidth',3); hold off;
    title(sprintf('gamma = %f\n', gammas(i)));
    drawnow;
end
plot(convflag)
    
%% Interesting! Let's make a cool plot
% if gamma < 0.55, then this doesn't seem to converge, if gamma > 0.55, it
% does.
% gamma = 0.55 + 0.001 - does not converge
% gamma = 0.55 + 0.003 - unclear
% gamma = 0.55 + 0.004 - unclear
% gamma = 0.55 + 0.005 - converges
niter = 1e6;
[x,hist,flag] = tpr.shifted(0.55,'maxiter',niter);
semilogy(hist); hold all; drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.003,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.004,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.005,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.01,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.1,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.2,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.55+0.3,'maxiter',niter);
semilogy(hist); drawnow;
%%
[x,hist,flag] = tpr.shifted(0.995,'maxiter',niter);
semilogy(hist); drawnow;

%% 
clf; hold all; set(gca,'YScale','log'); set(gca,'XScale','log');
gammas = fliplr([0.55 0.551 0.552 0.553 0.554 0.5545 0.555 0.56 0.65 0.995]);
ng = numel(gammas);
niter = 1e6;
nskip = 200;
for i=1:ng
    [x,hist,flag] = tpr.shifted(gammas(i),'maxiter',niter);
    hist = hist(1:nskip*floor(numel(hist)/nskip));
    %semilogy(hist); hold on;
    meanhist = mean(reshape(hist,nskip,numel(hist)/nskip));
    xmeanhist = nskip/2:nskip:numel(hist); % x coordinates of mean hist
    semilogy([1:(nskip/2-1) xmeanhist], [hist(1:(nskip/2-1))' meanhist]); 
    drawnow;
end

ylim([1e-5,1]);