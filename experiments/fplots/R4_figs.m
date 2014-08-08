function R4_figs
%% Load data
mats = load('../../tensors/mtm4.mat');
alpha = 0.99;
gamma = 1/2;
maxiter = 100;
n = 4;
v = ones(n,1)/n;
subset = {'R4_12','R4_17','R4_18','R4_19'};
%%
for i=1:numel(subset)
    %%
    %mat = mats.R4_mats{i};
    mat = subset{i};
    R = mats.(mat);
        
    tpr = tensorpr3(R,alpha,v);
    %smap = tetrahedron_plot_kappa(mat,tpr);
    %tetrahedron_plot_lambda2(mat,tpr);
        
    %[~,~,~,xhist1] = tpr.shifted(0);
    [~,~,~,xhistg] = tpr.shifted(gamma);
    [~,~,~,xhisti] = tpr.inverseiter;
    [~,~,~,xhistn] = tpr.newton;
    [~,~,~,xhistio] = tpr.innout;

    %%
    smap = tetrahedron_plot_lambda2(mat,tpr);
    myplot = @(x,varargin) plot3(x(1,:),x(2,:),x(3,:),varargin{:},'LineWidth',2.5);
    
        g.CameraPosition=[-3.2843 -10.4299 -6.2602];
          g.CameraTarget=[0.0990 -0.1826 0.2364];
        g.CameraUpVector=[0 0 1];
    g.CameraUpVectorMode='manual';
       g.CameraViewAngle=5.8808;
       set(gca,g);
       lighting gouraud
     camlight

    colorbar off;
    hold all;
    %%
    h = myplot(smap([v xhistg]),'g-');
    print(gcf,sprintf('%s-jac-shifted.eps',mat),'-depsc2');
    set(h,'Visible','off');
    %%
    h = myplot(smap([v xhistn]),'g-');
    print(gcf,sprintf('%s-jac-newton.eps',mat),'-depsc2');
    set(h,'Visible','off');
    %%
    h = myplot(smap([v xhistio]),'g-');
    print(gcf,sprintf('%s-jac-innout.eps',mat),'-depsc2');
    set(h,'Visible','off');
    %%
    h = myplot(smap([v xhisti]),'g-');
    print(gcf,sprintf('%s-jac-inverse.eps',mat),'-depsc2');
    set(h,'Visible','off');
    %%
    hold off;
    
end