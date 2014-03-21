%% Explore the counter-example in depth
% In this case, we use the simplex-plot ideas along with the
% counter-example to explore this case.

% These examples are for the CS 50th anniversary poster.,

alpha = 0.99;
gamma = 0.5;
v = ones(3,1)/3;

%% Requires a shift
R = [1/3  1/3  1/3  1/3  0  0.5  0  0  0
     1/3  1/3  1/3  1/3  0  0.5  1  1  1
     1/3  1/3  1/3  1/3  1  0  0  0  0];
 
tensorposter_study_figure('shift-required',R,alpha,v,gamma);
 
%% Grow in difference
R = [0.0000  0.0000  0.0000  1.0000  0.0000  0.5000  0.5000  1.0000  0.0000
     0.0000  0.0000  0.0000  0.0000  0.5000  0.5000  0.0000  0.0000  0.0000
     1.0000  1.0000  1.0000  0.0000  0.5000  0.0000  0.5000  0.0000  1.0000];
tensorposter_study_figure('step-growth',R,alpha,v,gamma); 

%% Jacobian large
R = [0.0000  0.0000  1/3  1/3  0.0000  0.0000  0.5000  0.5000  0.5000
     0.0000  0.0000  1/3  1/3  0.0000  0.0000  0.0000  0.0000  0.5000
     1.0000  1.0000  1/3  1/3  1.0000  1.0000  0.5000  0.5000  0.0000];
tensorposter_study_figure('large-jac',R,alpha,v,gamma);

