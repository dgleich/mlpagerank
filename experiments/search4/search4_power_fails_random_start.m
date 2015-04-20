
%%
f = @(T) ~get_outputs(@() tpr4_power(normout(reshape(T,3,27)')', ...
            0.99, [1/3;1/3;1/3]), 3); % get the flag output
search_tensors([3,3,3,3],f,'three4-did-not-converge-startnz', ...
    'filter',@(T) all(sum(reshape(T,3,27)) > 0),'startnz',27,'useperm',true);
    
