function [gamma,S] = li_gamma(P)
% LI_GAMMA Compute the gamma value from Li and Ng's paper
%
% [gamma,S] = li_gamma(P) directly computes the
% value of gamma and the associated set S for 
% a third order tensor. The tensor P can be given
% as either an n-times-n^2 array, or an 
% n-times-n-times-n tensor. 
%
% This computation works for larger tensors, but
% it must enumerate over all subsets, and so may
% be very expensive. 
%
% Example:
%   li_gamma([0 1 0 1; 1 0 1 0]) % should return 3
%  

if ismatrix(P)
    n = size(P,1);
    P = reshape(P,n,n,n);
end

if ndims(P) ~= 3
    error('li_gamma:invalidSize','gamma is only defined by 3x3x3 tensors');
end

n = size(P,1);

k = 1;
inds = {};
for i = 1:(n-1) 
    total = combnk(1:n, i);      % i-items subset 
    m = size(total, 1);
    
    for j = 1:m                  %j-th row contain the subset
        S = total(j, :);         %row vector
        Sbar = setdiff(1:n, S);  %row vector
        % first iterates all i3
        part1 = zeros(n,1);
        ind1 = 1;                %index for part1
        for i3 = 1:n            
            t1 = zeros(size(S,2),1);
            t2 = zeros(size(Sbar, 2),1);
            idx1 = 1;            %index for t1
            idx2 = 1;            %index for t2
            for i2 = S
                t1(idx1) = sum(P(Sbar,i2,i3));
                idx1 = idx1 + 1;
            end
            for i2 = Sbar
                t2(idx2) = sum(P(S,i2,i3));
                idx2 = idx2 + 1;
            end
            part1(ind1) = min(t1(1:(idx1-1))) + min(t2(1:(idx2-1)));
            ind1 = ind1 + 1;
        end
        
        % second iterates all i2
        part2 = zeros(n,1);
        ind2 = 1;               %index for part2
        for i2 = 1:n            
            t3 = zeros(size(S,2),1);
            t4 = zeros(size(Sbar,2),1);
            idx3 = 1;
            idx4 = 1;
            for i3 = S
                t3(idx3) = sum(P(Sbar,i2,i3));
                idx3 = idx3 + 1;
            end
            for i3 = Sbar
                t4(idx4) = sum(P(S,i2,i3));
                idx4 = idx4 + 1;
            end
            part2(ind2) = min(t3(1:(idx3-1))) + min(t4(1:(idx4-1)));
            ind2 = ind2 + 1;
        end
        tmp(k) = min(part1) + min(part2);
        inds(k) = {[i,j]};
        k = k + 1;
    end    
end

[gamma,ind] = min(tmp(1:(k-1)));
i = inds{ind}(1);
j = inds{ind}(2);

total = combnk(1:n, i);      % i-items subset 
S = total(j, :);         
