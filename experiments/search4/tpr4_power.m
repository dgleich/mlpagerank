function [x,hist,flag,ihist] = tpr4_power(R,alpha,v,varargin)

% SHIFTED Run the power method on a tensor PageRank problem
%
% x = shifted(P) solves with gamma=1/2, which may or may not
% converge.
%
% x = shifted(P,gamma) solves with a shift of gamma
% x = shifted(P,'key',value,'key',value,...)
% x = shifted(P,gamma,'key',value,'key',value,...)
%
% Valid additional parameters
%   'maxiter' : the maximum number of iterations, set to 1000
%   'tol' : the solution tolerance, set to 1e-8
%   'xtrue' : an optional true solution to get errors
%             to report for analytic descriptions

gamma=1;

p = inputParser;
p.addOptional('maxiter',1e4);
p.addOptional('tol',1e-8,@(x) isnumeric(x) && x<1 && x>0);
p.addOptional('xtrue',[]);
p.addOptional('randinit',false,@islogical);
p.addOptional('x0',[]);
p.parse(varargin{:});
opts = p.Results;
if nargout > 3, trackihist = 1; else trackihist = 0; end

% Extract data from obj
n = size(R,1);
a = alpha;

niter = opts.maxiter;
tol = opts.tol;

Gamma = 1 / (1+gamma);
xcur = zeros(n,1); % this allows us to keep v = 1/n :-)
xcur = xcur + v; 
if opts.randinit, xcur = rand(n,1); xcur=xcur/sum(xcur); end
if ~isempty(opts.x0), xcur = zeros(n,1) + opts.x0; end

hist = zeros(niter, 1);
if trackihist, ihist = zeros(n, niter); end

for i=1:niter
    % TODO make this iteration better
    y = a*(R*kron(kron(xcur, xcur), xcur)); 
    z = y * Gamma + Gamma*(1-sum(y))*v;
    xn = z + (1-sum(z))*xcur;

    if trackihist, ihist(:,i) = xn; end

    curdiff = norm(xcur - xn,1);
    curres = norm(a*(R*kron(kron(xn, xn), xn)) + (1-a)*v - xn, 1);
    hist(i) = curres;

    if ~isempty(opts.xtrue)
        hist(i) = norm(xn - opts.xtrue,inf);
    end

    % switch solutions
    xcur = xn;

    % check for termination
    if curres <= tol || curdiff <= tol/10;
        break;
    end                
end

hist = hist(1:i,:);
if trackihist, ihist = ihist(:, 1:i); end

if i == niter && curres > tol
    warning('tpr4_power:notConverged',...
        'did not converge after %i iterations to %e tolerance',...
        niter, tol);
    flag = 0;
else
    flag = 1;
end

x = xcur ./ sum(xcur);