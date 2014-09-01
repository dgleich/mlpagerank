M = combine_structs( ...
    load('../tensors/mtm3.mat'), ...
    load('../tensors/mtm4.mat'), ...
    load('../tensors/mtm6.mat'));

probs = [M.R3_mats M.R4_mats M.R6_mats];

n = 0;
i = 1;

fprintf('\\allowdisplaybreaks\n');

for pi=1:numel(probs)
    prob = probs{pi};
    R = M.(prob); % ick, but what to do?
    if size(R,1) ~= n
        if n > 0
            fprintf('\\end{align*} \n\n');
        end
        n = size(R,1);
        i = 1;
        fprintf('\\subsection{%i $\\times$ %i $\\times$ %i}\n', n,n,n);
        fprintf('\\begin{align*} \n');
    end
    
    fprintf('\\mR_{%i,%i} & = \\left[ \\scriptstyle', n, i);
    print_tensor3(spones(R));
    fprintf('\\right] \\\\ \n');
    
    i = i+1;
end
fprintf('\\end{align*} \n\n');
    