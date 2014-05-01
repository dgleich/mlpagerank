%% Try and trace down this convergence result
a = 1/2-0.001;
%a = 1/2 - 0.3;
f = (1-a);
fk = [];
for i=1:25
    theta = a/((1-2*a)^2);
    f = a*f^2/((1-2*a)^2);
    fk(i) = f;
    s = '';
    %gk = (1-a)^2*(2*a)^(2^i);
    
    %c = (a*theta - 1)/(2*theta);
    c = max(1-2*a,2*a);
    gk = min(theta*(c)^(2^i),(1/4)^(i-1)*(1-a)^2*a);
    if gk < f, s = '****'; end
    fprintf('%4i  %.20f  %.20f  %.20f %s\n', i, f, gk, (1/4)^(i-1)*(1-a)^2*a, s);
    %fprintf('%4i  %.20f  %.20f\n', i, f, (1-2*a)^2/(a)*((1-2*a)^(2^(i))));
end

% Try and determine the rate of convergence
rates = fk(2:end) ./ fk(1:end-1).^2;
rates = rates(isfinite(rates));
fprintf('rate = %.6f   a = %.6f   guess = %.6f\n', rates(end), a, theta);
