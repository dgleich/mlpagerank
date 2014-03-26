%% Load data
mats = load('../../tensors/mtm3.mat');

%%
for i=1:numel(mats.R3_mats)
    mat = mats.R3_mats{i};
    R = mats.(mat);
    %triangle_plot_lambda2(mat,R,0.99,ones(3,1)/3,1/2);
    triangle_plot_kappa(mat,R,0.99,ones(3,1)/3,1/2);
    mat
    pause;
end