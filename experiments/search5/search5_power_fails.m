f = @(T) ~get_outputs(@() tpr5_power(normout(reshape(T,2,16)')', ...
            0.99, [1/2;1/2]), 3); % get the flag output
search_tensors([2,2,2,2,2],f,'two5-did-not-converge', ...
    'filter',@(T) all(sum(reshape(T,2,16)) > 0),'useperm',true);

%% This completed and didn't find any examples

