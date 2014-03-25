function genSymScript(R)
% This function generates a matlab script file called 'mySym.m'
% The files contains symbolic equation of the tensorPageRank problem.
fid = fopen('mySym.m', 'w');
n = min(size(R));
fprintf(fid, 'function Y = mySym(R)\n');
s = 'syms ';
t = '';
for i = 1:n
    t = [t, 'x', int2str(i), ' '];
end
s = [s, t, '\n'];
fprintf(fid, s);
fprintf(fid, 'alpha = 0.99;\n');
s = ['v = 1/', int2str(n), ';\n'];
fprintf(fid, s);
sol = 'sol_';
t = '[';
for i = 1:n-1
    t = [t, sol, 'x', int2str(i), ', '];
end
t = [t, sol, 'x', int2str(n), ']', ' = ...\n'];
fprintf(fid, t);
fprintf(fid, '\tsolve(');
for i = 1:n
    tmp = ['x', int2str(i), ' == ', 'alpha*(0'];
    for j = 1:n^2
        if (R(i, j) ~= 0)
            jindex = mod(j, n);
            if (jindex ~= 0)
                iindex = (j - jindex) / n + 1;
            else
                jindex = n;
                iindex = j / n;
            end
            tmp = [tmp, ' + ', 'R(', int2str(i), ',', int2str(j), ')*', 'x', ...
                int2str(iindex), '*x', int2str(jindex)]; % map j to x index
        end
    end
    if (i ~= n)
        tmp = [tmp, ') + (1-alpha)*v, ...\n\t\t  '];
        fprintf(fid, tmp);
    else
        tmp = [tmp, ') + (1-alpha)*v);\n'];
        fprintf(fid, tmp);
    end
end
s = ['Y = double(['];
for i = 1:n
    s = [s, sol, 'x', int2str(i), ' '];
end
s = [s, ']'');\n'];
fprintf(fid, s);
% post-processing
p = ['k = 1;\n'];
p = [p, 'X = zeros(size(Y));\n'];
p = [p, 'for j = 1:size(Y, 2)\n\t'];
p = [p, 'if (min(real(Y(:, j))) >= 0 && norm(max(imag(Y(:, j)))) <= eps(1) && norm(sum(Y(:, j))-1) <= eps(1))\n\t\t'];
p = [p, 'X(:, k) = Y(:, j);\n\t\t'];
p = [p, 'k = k + 1;\n\t'];
p = [p, 'end\nend\n'];
p = [p, 'Y = X(:, 1:k-1);\n'];
fprintf(fid, p);
fprintf(fid, 'end\n');
fclose(fid);
end
