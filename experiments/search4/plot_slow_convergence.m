Ts = load('three4-did-not-converge-found');

fn = fieldnames(Ts);
fi = 32;
Tname = fn{fi};
T = Ts.(Tname);

[~,hist,flag,xhist] = tpr4_power(normout(reshape(T,3,27)')', ...
                0.99, [1/3;1/3;1/3], 'maxiter', 1e5);

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
set_figure_size([5,3.5]);
box off;
axes('position',[0.65 0.6 0.25 0.25]);
reshist = zeros(size(xhist,2),1);

semilogy(hist,'k.-');
%xlim([1,5]);
ylim([1e-10,1]);
xlim([1,size(xhist,2)])
xlabel('Iteration');
ylabel('Residual');
box off;

srate = floor(size(xhist,2)/2); 
rates = exp(diff(log(hist)));
frate = mean(rates(end-srate:end));
text(srate,hist(srate),sprintf('Rate\n%.4f',frate), 'VerticalAlignment','bottom');

print('three4-slow-converge.eps','-depsc2');

% !cp three4-slow-converge.eps ~/Dropbox/publications/tensorpr-shared/figures/

%% Write out the entries
nzs = find(T);
for nzi=1:numel(nzs)
    [i1,i2,i3,i4] = ind2sub(size(T),nzs(nzi));
    fprintf('\\cA(%i,%i,%i,%i) = 1 \\\\\n',i1,i2,i3,i4);
end