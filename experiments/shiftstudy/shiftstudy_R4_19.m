%% Check if R4_1 converges for any shift.

load ../../tensors/mtm4

R = R4_19;
n = size(R,1);
alpha = 0.99;
tpr = tensorpr3(R,alpha);
ng = 100;
gammas = linspace(0,1,ng+2);
gammas = gammas(2:end-1); % remove end points
convflag = zeros(ng,1);

%%
for i=1:ng
    [x,hist,flag] = tpr.shifted(gammas(i),'maxiter',1000);
    convflag(i) = 1;
end
plot(convflag)
    
%% Nope!
%% Try a very long time

xg=tpr.shifted('maxiter',100000);
%%
xg=tpr.shifted(2/3,'maxiter',100000);