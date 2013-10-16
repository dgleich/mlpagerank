addpath('../../../');
alpha = 0.99;
v = ones(3,1)./3;
R1 = [1/3 1/3 1/3 1/3 0 0 0 0 0;...
      1/3 1/3 1/3 1/3 0 1/2 1 0 1;...
      1/3 1/3 1/3 1/3 1 1/2 0 1 0];
R3 = [0 0 0 1 0 1/2 1/2 1 0;...
      0 0 0 0 1/2 1/2 0 0 0;...
      1 1 1 0 1/2 0 1/2 0 1];
%% fig 4.3(a)
name = 'NoNeedShift';  
triangle_plot_kappa(name, R3, alpha, v, 1/2);  
%% fig 4.3(b)
name = 'NoNeedShift';
triangle_plot_lambda2(name, R3, alpha, v, 1/2);
%% fig 4.3(c)
name = 'NeedShift';
triangle_plot_kappa(name, R1, alpha, v, 1/2);  
%% fig 4.3(d)
name = 'NeedShift';
triangle_plot_lambda2(name, R1, alpha, v, 1/2);