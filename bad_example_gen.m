function R = bad_example_gen(n)
% R = bad_example_gen(n)
% this function takes the dimension n as the argument 
% and returns a nxn^2 matrix R
Rs = {};
for i = 1:n
    Rs{i} = full(spones(sprand(n, n, log(i)/n)));
    Rs{i}(:, n) = Rs{i}(:, n) + (sum(Rs{i}, 2) == 0);
    A = Rs{i}';
    % normout ith cell matrix
    for j = 1:n
        A(:, j) = A(:, j) ./ sum(A(:, j));
    end
    Rs{i} = A;
end
R = [Rs{:}];
end