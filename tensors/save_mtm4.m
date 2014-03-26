clear
load old/4x4x4_not_converge_non_shift
R4_mats = {};
for i=1:18
    mat = sprintf('R4_%i',i);
    R = eval(sprintf('R%i',i));
    assignin('base',mat,R);
    R4_mats{end+1} = mat;
end
%%
for i = 1:numel(R4_mats) 
    mat = R4_mats{i};
    R = eval(mat);
    R4_Properties.(mat).gamma = li_gamma(R);
    R4_Properties.(mat).alpha99.gamma = li_gamma(tensorpr3(R,0.99).markov());
    R4_Properties.(mat).alpha95.gamma = li_gamma(tensorpr3(R,0.95).markov());
    R4_Properties.(mat).alpha90.gamma = li_gamma(tensorpr3(R,0.90).markov());
    R4_Properties.(mat).alpha85.gamma = li_gamma(tensorpr3(R,0.85).markov());
    R4_Properties.(mat).alpha70.gamma = li_gamma(tensorpr3(R,0.70).markov());
end
%%
save('mtm4.mat',R4_mats{:},'R4_mats','R4_Properties'); 