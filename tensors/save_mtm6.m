clear
z = 1;
load old/6x6x6_not_converge_non_shift
R6_mats = {};
for i=1:2
    mat = sprintf('R6_%i',z); z = z+1;
    R = eval(sprintf('R%i',i));
    assignin('base',mat,R);
    R6_mats{end+1} = mat;
end

load old/6x6x6_slow_converge_shift
for i=1:2
    mat = sprintf('R6_%i',z); z = z+1;
    R = eval(sprintf('T%i',i));
    assignin('base',mat,R);
    R6_mats{end+1} = mat;
end

load old/counterExample
mat = sprintf('R6_%i',z); z = z+1;
R = eval('T');
assignin('base',mat,R);
R6_mats{end+1} = mat;

%%
for i = 1:numel(R6_mats) 
    mat = R6_mats{i};
    R = eval(mat);
    R6_Properties.(mat).gamma = li_gamma(R);
    R6_Properties.(mat).alpha99.gamma = li_gamma(tensorpr3(R,0.99).markov());
    R6_Properties.(mat).alpha95.gamma = li_gamma(tensorpr3(R,0.95).markov());
    R6_Properties.(mat).alpha90.gamma = li_gamma(tensorpr3(R,0.90).markov());
    R6_Properties.(mat).alpha85.gamma = li_gamma(tensorpr3(R,0.85).markov());
    R6_Properties.(mat).alpha70.gamma = li_gamma(tensorpr3(R,0.70).markov());
end
%%
save('mtm6.mat',R6_mats{:},'R6_mats','R6_Properties'); 

%% Add exact solutions
add_exact('mtm6.mat');