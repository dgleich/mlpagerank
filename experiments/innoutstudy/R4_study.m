load ../../tensors/mtm4.mat
% This is the set of problems from R4 that did not converge
% but some of them may ...
innout_set = {'R4_3','R4_12','R4_17'};
alpha = 0.99;
for i=1:numel(innout_set)
    mat = innout_set{i};
    R = eval(mat);
    tpr = tensorpr3(R,alpha);
    x = tpr.innout('maxiter',2000);
    mat
    norm(tpr.residual(x))
end
%%
% R4_17 seems to be the hardest
% Except shifted works!
% It has multiple solutions.