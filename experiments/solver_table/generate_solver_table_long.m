%% compute a table of the performance for all the solvers are the test problems

M = combine_structs( ...
    load('../../tensors/mtm3.mat'), ...
    load('../../tensors/mtm4.mat'), ...
    load('../../tensors/mtm6.mat'));

probs = [M.R3_mats M.R4_mats M.R6_mats];
Props = combine_structs(M.R3_Properties, M.R4_Properties, M.R6_Properties);
%%
alphas = [0.7 0.85 0.9 0.95 0.99];
methods = {'power','shifted','innout','inverseiter','newton'};
niters = [100000,100000,10000,10000,10000];

%% Solver table 1, number of problems that converge with default settings
fprintf('%10s  %5s  %15s  %10s\n', 'problem', 'alpha', 'method','result');
warning('off','tensorpr3:notConverged');
results = zeros(numel(probs),numel(alphas),numel(methods));
solutions = cell(numel(probs),numel(alphas));
for pi=1:numel(probs)
    prob = probs{pi};
    R = M.(prob); % ick, but what to do?
    for ai=1:numel(alphas)
        tpr = tensorpr3(R,alphas(ai));
        X = zeros(size(R,1),numel(methods));
        for mi=1:numel(methods)
            [X(:,mi),~,flag] = tpr.(methods{mi})('maxiter',niters(mi));
            if flag==1, result='success'; else result='failed ***'; end;
            fprintf('%10s  %5.2f  %15s  %10s\n',...
                prob, alphas(ai), methods{mi}, result);
            results(pi,ai,mi) = flag;
        end
        solutions{pi,ai} = X;
    end
end
warning('on','tensorpr3:notConverged');

save 'solver_long.mat' M probs Props alphas methods solutions results
%%
load solver_long

%% Prepare a table of the results
% $\alpha$ & $n$ & Method & \\
%          &     & Power & Shifted & Newton & Innout & Inverse \\
%  0.70    & 3   & 
for ai = 1:numel(alphas)
    probsize = zeros(numel(probs),1);
    
    for pi = 1:numel(probs)
        probsize(pi) = size(solutions{pi,ai},1);
    end
    probresults = squeeze(results(:,ai,:)); % this is a nproblems-by-nmethods matrix
    ns = sort(unique(probsize));
    
    % aggregate results by
    for si = 1:numel(ns)
        n = ns(si);
        nprobs = sum(probsize==n);
        
        if si==1
            fprintf(' %.2f ', alphas(ai));
        else
            fprintf('      ');
        end
        fprintf(' & ');
        fprintf(' %i ', n);
        
        for mi=1:numel(methods)
            fprintf(' & %i ', sum(probresults(probsize==n,mi)));
        end
        
        fprintf(' \\\\ \n');
    end
    
    % write out the summary
    
    fprintf('       &    ');
    
    for mi=1:numel(methods);
        fprintf(' & %i ', sum(probresults(:,mi)));
    end
    
    fprintf(' \\\\ \n');
end
        