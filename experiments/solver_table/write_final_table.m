L = load('solver_long');
D = load('solver_defaults');
%%
assert(all(L.alphas == D.alphas));
assert(all(strcmp(L.methods,D.methods)));
assert(all(strcmp(L.probs,D.probs)));

%%
load 'solver_defaults'; % This is the first set of columns

%% Prepare a table of the results
% $\alpha$ & $n$ & Method (Default) & ... & Method (Long) \\
%          &     & Power & Shifted & Newton & Innout & Inverse & 
%                      Power & Shifted & Newton & Innout & Inverse \\
%  0.70    & 3   & 
for ai = 1:numel(alphas)
    probsize = zeros(numel(probs),1);
    
    for pi = 1:numel(probs)
        probsize(pi) = size(solutions{pi,ai},1);
    end
    probresults = squeeze(results(:,ai,:)); % this is a nproblems-by-nmethods matrix
    Lprobresults = squeeze(L.results(:,ai,:)); % this is a nproblems-by-nmethods matrix
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
        
        for mi=1:numel(methods)
            fprintf(' & %i ', sum(Lprobresults(probsize==n,mi)));
        end
        
        if si==numel(ns)
            % last one
            fprintf(' \\\\ \\addlinespace \n');
        else
            fprintf(' \\\\ \n');
        end
    end
    
    % write out the summary
    
    fprintf('       &    ');
    
    for mi=1:numel(methods);
        fprintf(' & %i ', sum(probresults(:,mi)));
    end
    
    for mi=1:numel(methods);
        fprintf(' & %i ', sum(Lprobresults(:,mi)));
    end
    
    fprintf(' \\\\ \\midrule \n');
end
        