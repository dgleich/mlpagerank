function [Y,Yall,eqs,sol] = solve_tensorpr3_symbol(R,astr)
% SOLVE_TENSORPR3_SYMBOL Solve a tensorpr3 problem symbolically.
%
% [Y,Yall] = solve_tensorpr3_symbol(R,astr) uses the Matlab symbolic
% toolbox to compute all solutions to the polynomial system of equations: 
%
%   astr*R*kron(x,x) + (1-astr)/n = x
%
% where astr is a string that evaluates alpha, and the entire expression is
% written symbolically. The matrix Yall returns all solutions computed,
% whereas Y only returns those that satisfy the tensorpr3 equation -- that
% is, they have sum=1, and are non-negative.
%
% [Y,Yall,eqs,sol] = solve_tensorpr3_symbol(R,astr) also returns the vector of
% equations along with the symbolic solution structure. This is useful for
% degenerate problems.

%%
n = size(R,1);
strs = cell(n,1);
for j=1:n
    curstr = cell(n^2+2,1);
    s = 1;
    for k1=1:n
        for k2=1:n
            aij = rats(R(j,s));
            curstr{s} = sprintf('(%s)*(%s)*%s*%s + ',astr,aij,['x' k1+'0'],['x' k2+'0']);
            s = s+1;
        end
    end
    curstr{n^2+1} = sprintf(' (1-%s)/%i = ',astr,n);
    curstr{n^2+2} = ['x' j+'0'];
    strs{j} = char([curstr{:}]);
end
eqs = strs;
vars = cell(n,1);
for i=1:n
    vars{i} = sprintf('x%i',i);
end
sol = solve(strs{:},vars{:});
solcell = struct2cell(sol);
if ~all(cellfun(@(x) isempty(symvar(x)),solcell))
    % that means the problem is denengerate
    fprintf('Degenerate problem detected, pruning solution spaces\n');
    degenerate = zeros(size(solcell{1},1),1);
    for j=1:numel(solcell)
        for i=1:numel(solcell{j})
            degenerate(i) = degenerate(i) | ~isempty(symvar(solcell{j}(i)));
        end
    end
    for j=1:numel(solcell)
        solcell{j} = solcell{j}(~degenerate);
    end
end
    
Y = cell2mat(cellfun(@double,solcell,'UniformOutput',false)')';
Yall = Y;
Y = Y(:,~any(abs(imag(Y)) >= sqrt(eps)));
Y = real(Y);
Y = Y(:,~any(Y <0));
Y = Y(:,~any(Y >1));
Y = Y(:,sum(Y) >= 1-sqrt(eps) & sum(Y) <= 1 + sqrt(eps)); 