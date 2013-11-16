function [gamma,S] = computeGamma(R, alpha, v)
% [gamma,S] = computeGamma(R, alpha, v)
% argument R is the original transition tensor
% alpha is the damping factor
% gamma will be intractable for large dimensions
% the order of the tensor should be 3
n = size(R, 1); % n is the dimension
if nargin < 3
    v = ones(n, 1)./n;
end
e = ones(n,1);
T = alpha .* R + (1-alpha).*v *kron(e',e');
P = zeros(n,n,n);
tmp = zeros(2^n - 2, 1);% totally (2^n -2) non-empty subset
k = 1;
for i = 1:n
    P(:,i,:) = T(:, (i-1)*n+1:i*n);
end
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

end