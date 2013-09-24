% I want to test a few properties on reductions of the two-step chain
% inspired by the Quadratic DS paper by Y. Rabinowich.
for i=1:1000
    n = 5;
    y = rand(n^2,1); y = y/sum(y); % generate an n^2 vector
    x = rand(n,1); x = x/sum(x);
    sy = sum(reshape(y,5,5),2);
    if norm(kron(x,x) - kron(sy,sy),1) >= norm(kron(x,x) - y,1)
        fprintf('found 1\n');
        norm(kron(x,x) - y,1)
        norm(kron(x,x) - kron(sy,sy),1)
        break;
    end
end
% Nothing here worked, damn.