cmap = hot(32);
cmap = flipud(cmap); % make white first
cmap = cmap(9:24,:); % rescale until red

% cmap2 = copper(48);
% cmap2 = cmap2(1:32,:);
% cmap2(:,1) = 0;
% cmap2(:,2) = cmap2(:,2)*0.5;
% cmap2(:,3) = cmap2(:,3)*0.8;

cmap2 = gray(32);
cmap2 = cmap2(1:32,:);
cmap2(:,1) = 0.1;
cmap2(:,2) = cmap2(:,2)*0.6;
cmap2(:,3) = cmap2(:,3)*0.9;


cmap = [cmap2; cmap];
colormap(cmap);
caxis([0 1.5]);