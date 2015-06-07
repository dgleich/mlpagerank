function [x,hist] = mlpr_large(T,alpha,varargin)
% MLPR_LARGE Solve a large MLPR problem given 

parser = inputParser;
parser.addOptional('maxiter',1000);
parser.addOptional('tol',1e-8);
parser.addOptional('v',[]);
parser.addOptional('shift',0);
parser.parse(varargin{:});
opts = parser.Results;

n = max(max(T));

colnorms = norms(T,n);

if issempty(opts.v) 
    v = ones(n,1)/n;
else
    v = opts.v;
end

maxiter=opts.maxiter;
tol=opts.tol;

x = full(v);
hist = zeros(maxiter,1);

for iter=1:maxiter
    y = mult(n,T,norms,x);
    y = y + (1-sum(y))*v;
    y = alpha*y + (1-alpha)*v;
    
    hist(iter) = norm(y-x,1);
    if hist(iter) < tol
        break
    end
    
    x = (1-shift)*y + shift*x;
end

function colnorms=norms(T,n)

S = sparse(T(:,2),T(:,3),1,n,n);
norms = zeros(size(T,1),1);

inds = sub2ind(T(:,2),T(:,3),n,n)
colnorms = S(inds);

function y=mult(n,T,norms,x)
y = zeros(n,1);
for i=1:size(T,1)
    ci = T(i,1);
    cj = T(i,2);
    ck = T(i,3);
    denom = norms(i);
    y(ci) = y(ci) + x(cj)*x(ck)/denom;
end

