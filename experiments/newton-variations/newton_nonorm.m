function [x, hist, flag, ihist] = newton_ls(obj, varargin)
    % NEWTON Solve the tensorpr3 iteration using Newton's method

    p = inputParser;
    p.addOptional('maxiter',1e3);
    p.addOptional('tol',1e-8,@(x) isnumeric(x) && x<1 && x>0);
    p.addOptional('xtrue',[]);
    p.addOptional('x0',[]);
    p.addOptional('randinit',false,@islogical);
    p.parse(varargin{:});
    opts = p.Results;
    if nargout > 3, trackihist = 1; else trackihist = 0; end

    % Extract data from obj
    R = obj.R;
    n = size(R,1);
    a = obj.alpha;
    v = obj.v;

    niter = opts.maxiter;
    tol = opts.tol;
    %xcur = zeros(n,1) + v;
    xcur = zeros(n,1) + (1-a)*v;
    if opts.randinit, xcur = rand(n,1); xcur=xcur/sum(xcur); end
    if ~isempty(opts.x0), xcur = zeros(n,1) + opts.x0; end

    hist = zeros(niter, 1);
    if trackihist, ihist = zeros(n, niter); end

    I = eye(n);
    for i = 1:niter
        % This iteration is equivalent to Newton's method
        % for a fixed point g(x) - x = 0.
        % where the iteration is (J(x) - I) x_{k+1} = J(x)*x -
        % g(x),
        % which for this form


        A = a*R*(kron(xcur, I) + kron(I, xcur)) - I;
        b = a*R*kron(xcur, xcur) - (1-a)*v; 

        xn = A\b;
                
        if trackihist, ihist(:,i) = xn; end

        curdiff = norm(xcur - xn,1) + (1-sum(xn));
        curres = norm(obj.residual(xn), 1) + (1-sum(xn));
        hist(i) = curres;

        if ~isempty(opts.xtrue), hist(i) = norm(xn - opts.xtrue,inf); end

        xcur = xn;

        if curres <= tol || curdiff <= tol/10;
            break
        end 
    end

    hist = hist(1:i, :);
    if trackihist, ihist = ihist(:, 1:i); end
    if i == niter && curres > tol
        warning('tensorpr3:notConverged',...
            'did not converge after %i iterations to %e tolerance',...
            niter, tol);
        flag = 0;
    else
        flag = 1;
    end

    x = xcur;
end                

