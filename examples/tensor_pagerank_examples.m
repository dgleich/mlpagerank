%% Example for paper

load tensors/example1
n = size(R,1);

S = tensorpr3(R,0.85);
x = S.solve('tol',1e-16);
x
