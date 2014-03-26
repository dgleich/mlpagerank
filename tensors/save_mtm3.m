%% Save the initial set of tensors into mtm3
R3_mats = {};
R3_1 = [1/3 1/3 1/3  1/3 0  0   0 0 0
      1/3 1/3 1/3  1/3 0 1/2  1 0 1
      1/3 1/3 1/3  1/3 1 1/2  0 1 0];
R3_mats{end+1} = 'R3_1'; 

R3_2 = [0 0 0  1  0  1/2  1/2 1 0
        0 0 0  0 1/2 1/2   0  0 0
        1 1 1  0 1/2  0   1/2 0 1];  
R3_mats{end+1} = 'R3_2'; 

R3_3 = [ 0 1/2 0  1/2  0  1    1/2 1/2  0   
         0  0  0   0  1/2 0     0  1/2  0   
         1 1/2 1  1/2 1/2 0    1/2  0   1  ];
R3_mats{end+1} = 'R3_3'; 

R3_4 = [0.0000  0.0000  1/3  1/3  0.0000  0.0000  0.5000  0.5000  0.5000
        0.0000  0.0000  1/3  1/3  0.0000  0.0000  0.0000  0.0000  0.5000
        1.0000  1.0000  1/3  1/3  1.0000  1.0000  0.5000  0.5000  0.0000];
R3_mats{end+1} = 'R3_4';    

%% Compute the 
for i = 1:numel(R3_mats) 
    mat = R3_mats{i};
    R = eval(mat);
    R3_Properties.(mat).gamma = li_gamma(R);
    R3_Properties.(mat).alpha99.gamma = li_gamma(tensorpr3(R,0.99).markov());
    R3_Properties.(mat).alpha95.gamma = li_gamma(tensorpr3(R,0.95).markov());
    R3_Properties.(mat).alpha90.gamma = li_gamma(tensorpr3(R,0.90).markov());
    R3_Properties.(mat).alpha85.gamma = li_gamma(tensorpr3(R,0.85).markov());
    R3_Properties.(mat).alpha70.gamma = li_gamma(tensorpr3(R,0.70).markov());
end
%% Save
save('mtm3.mat',R3_mats{:},'R3_mats','R3_Properties'); 