Ts = load('three4-did-not-converge-found');

Tfail = {};
fn = fieldnames(Ts);
for fi=1:numel(fn)
    Tname = fn{fi};
    T = Ts.(Tname);
    [~,~,flag] = tpr4_power(normout(reshape(T,3,27)')', ...
                    0.99, [1/3;1/3;1/3], 'maxiter', 1e5);
    if ~flag
        Tfail{end+1} = Tname;
    else
        fprintf('%4s passed\n', Tname);
    end
end

