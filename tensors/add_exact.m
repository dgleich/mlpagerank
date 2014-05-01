%% add the exact solutions to a file
function add_exact(filename)

%%
%filename = 'mtm3.mat';
D = load(filename);

%%
vars = fieldnames(D);
matsvar = vars(~cellfun(@isempty, strfind(vars,'_mats')));
mats = D.(matsvar{1});
prefix = strtok(matsvar,'_');
prefix = prefix{1};
%%
alphas={'70','85','90','95','99','100'};
for i=1:numel(mats)
    mat = D.(mats{i});
    % create the strings for each equation
    for ai=1:numel(alphas)
        alpha = alphas{ai};
        astr = sprintf('%s/100',alpha);
        [Y,Yall] = solve_tensorpr3_symbol(mat,astr);
        fprintf('Problem %10s with a=%10s had %2i stochastic solutions of %2i total sol\n',...
            mats{i}, astr, size(Y,2), size(Yall,2));
        if strcmp(alpha,'100')
            D.(sprintf('%s_Properties',prefix)).(mats{i}).sols = Y;
            D.(sprintf('%s_Properties',prefix)).(mats{i}).all_sols = Yall;
        else
            D.(sprintf('%s_Properties',prefix)).(mats{i}).(sprintf('alpha%s',alpha)).sols = Y;
            D.(sprintf('%s_Properties',prefix)).(mats{i}).(sprintf('alpha%s',alpha)).all_sols = Yall;
        end
    end
end
%%
save(filename,'-struct','D');

