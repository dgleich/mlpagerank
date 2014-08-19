function print_tensor3(R)
R = full(R);
n = size(R,1);
fprintf('\\begin{array}');
fprintf('{');
pattern = repmat('c',1,n);
fprintf('%s',pattern);
for i=1:n-1
    fprintf('|%s',pattern);
end
fprintf('}\n');
for i=1:n
    for j=1:n^2-1
        fprintf('%i & ', R(i,j));
    end
    fprintf('%i \\\\\n', R(i,end));
end
fprintf('\\end{array}\n');