clc;
fprintf('Start 6x6x6 bad examples...\n');
n = 6;
N = 5e6; % test millions of bad-example candidates
fid1 = fopen('bad_example_6x6x6_shift.txt', 'a');
fid2 = fopen('bad_example_6x6x6_new.txt', 'a');
for i = 1:N
    % generate bad_example candidates
    R = bad_example_gen(n);
    alpha = 0.99;
    e = ones(n, 1);
    v = e ./ n;
    % shifted method
    [~, ~, flag1, ~, ~, ~, ~] = tensorRank(alpha, R, v, v);
    % newUpdate method
    [~, ~, ~, flag2] = newUpdate(alpha, R, v, v);
    if flag1 > 0
        dlmwrite('bad_example_6x6x6_shift.txt,', R, '-append', 'delimiter',...
            '\t', 'precision', '%.4f');
        fprintf(fid1, '\n');
    end
    if flag2 > 0
        dlmwrite('bad_example_6x6x6_new.txt', R, '-append', 'delimiter',...
            '\t', 'precision', '%.4f');
        fprintf(fid2, '\n');
    end
end
fclose(fid1);
fclose(fid2);
fprintf('Done!\n');
