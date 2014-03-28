cmap = hot(32);
cmap = flipud(cmap); % make white first
cmap = cmap(9:24,:); % rescale until red

cmap = [cmap];
colormap(cmap);
caxis([0 1]);