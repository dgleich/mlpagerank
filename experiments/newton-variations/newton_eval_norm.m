%% Evaluate whether or not we should normalize newton's method
% In these problems, we just check if we should use Newton's method
% with or without explicit stochastic normalization
M = combine_structs( ...
    load('../../tensors/mtm3.mat'), ...
    load('../../tensors/mtm4.mat'), ...
    load('../../tensors/mtm6.mat'));

probs = [M.R3_mats M.R4_mats M.R6_mats];
Props = combine_structs(M.R3_Properties, M.R4_Properties, M.R6_Properties);

alphas = [0.2 0.4 0.499 0.7 0.85 0.9 0.95 0.99];
methods = { {'snewton', @(tpr) newton(tpr)}, 
            {'nonorm', @(tpr) newton_nonorm(tpr)},
            {'nonorm0', @(tpr) newton_nonorm(tpr,'x0',0)} };

%%
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
            mname = methods{mi}{1};
            mfunc = methods{mi}{2};
            [X(:,mi),~,flag] = mfunc(tpr);
            if flag==1, result='success'; else result='failed ***'; end;
            fprintf('%10s  %5.2f  %15s  %8.1e  %8.1e  %-10s\n',...
                prob, alphas(ai), mname, norm(tpr.residual(X(:,mi)),1), sum(X(:,mi)), result);
            results(pi,ai,mi) = flag;
        end
        solutions{pi,ai} = X;
    end
end

%%
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

%% Conclusion
% Newton without normalization will not solve any problem if alpha > 0.5
% So we shouldn't use this one.
