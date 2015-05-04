f = @(T) ~get_outputs(@() tpr4_power(normout(reshape(T,2,8)')', ...
            0.99, [1/2;1/2]), 3); % get the flag output
search_tensors([2,2,2,2],f,'two4-did-not-converge');

%%
f = @(T) ~get_outputs(@() tpr4_power(normout(reshape(T,3,27)')', ...
            0.99, [1/3;1/3;1/3]), 3); % get the flag output
search_tensors([3,3,3,3],f,'three4-did-not-converge');

