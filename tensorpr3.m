classdef tensorpr3
    properties
        R
        alpha
        v
    end
    
    methods
        function obj = tensorpr3(R, alpha, v)
            % TENSORPR3 Create a TensorPR3 problem
            % A tensor PR3 problem solves:
            %   x = alpha*R*kron(x,x) + (1-alpha) v
            % The matrix R is required.
            % The value alpha defaults to 1/2
            % and the value v defaults to 1/n, where R is n-by-n^2
            % The matrix R must be column stochastic.
            
            n = size(R,1);
            if nargin < 2
                alpha = 1/2;
            end
            if nargin < 3
                v = ones(n, 1) ./ n;
            end
            % error checking
            if norm(ones(1,n)*v - 1, 1) > eps(1)
                error('input vector v is not stochastic.');
            else
                v = v ./ sum(v);
            end
            
            if min(ones(1, n) * R) ~= 1
                error('input matrix R is not column stochastic.');
            end
            
            obj.R = R;
            obj.alpha = alpha;
            obj.v = v;
        end
        
        function J = jacobian(obj,x,gamma)
            % JACOBIAN return the Jacobian of the problem at x with shift
            % gamma
            
            n = size(obj.R,1);
            I = eye(n);
            J = obj.alpha*gamma*obj.R*(kron(x,I) + kron(I,x)) + (1-gamma)*I;
        end
        
        function r = residual(obj,x)
            r = obj.alpha * (obj.R * kron(x, x)) + (1-obj.alpha) * obj.v - x;
        end
        
        function [x,hist,flag] = solve(obj,varargin)
            % SOLVE Run the power method on a tensor PageRank problem
            %
            % x = solve(P) solves with gamma=1/2, which may or may not
            % converge.
            %
            % x = solve(P,gamma) solves with a shift of gamma
            % x = solve(P,'key',value,'key',value,...)
            % x = solve(P,gamma,'key',value,'key',value,...)
            %
            % Valid additional parameters
            %   'maxiter' :
            %   'tol' :
            %   'xtrue' :
            % TODO, finish the description
            
            
            if nargin>1 && isnumeric(varargin{1})
                gamma = varargin{1};
                varargin = varargin(2:end);
            else
                gamma = 0.5;
            end
            
            p = inputParser;
            p.addOptional('maxiter',1e5);
            p.addOptional('tol',1e-10);
            p.addOptional('xtrue',[]);
            p.parse(varargin{:});
            opts = p.Results;
            
            % Extract data from obj
            R = obj.R;
            n = size(R,1);
            a = obj.alpha;
            v = obj.v;
            
            niter = opts.maxiter;
            tol = opts.tol;
            
            Gamma = 1 / (1+gamma);
            xcur = zeros(n,1); % this allows us to keep v = 1/n :-)
            xcur = xcur + v; 
            
            hist = zeros(niter, 1);
            
            for i=1:niter
                % TODO make this iteration better
                y = a*(R*kron(xcur, xcur)); % make sure Matlab does it 
                z = y * Gamma + Gamma*(1-sum(y))*v;
                xn = z + (1-sum(z))*xcur;
                
                curdiff = norm(xn - xcur, 1);
                hist(i) = curdiff;
                if ~isempty(opts.xtrue)
                    hist(i) = norm(xn - opts.xtrue,inf);
                end
                
                % check for termination
                if curdiff <= tol
                    break;
                end
                
                % switch solutions
                xcur = xn;
            end
            
            hist = hist(1:i,:);
            if i == niter && curdiff > tol
                warning('did not converge');
                flag = 0;
            else
                flag = 1;
            end
            
            x = xn / sum(xn);
        end
        
        function [x, hist, flag] = solven(obj, varargin)
            % Solven solve the tensorpr3 iteration using non-shift method
            
            p = inputParser;
            p.addOptional('maxiter',1e5);
            p.addOptional('tol',1e-8);
            p.addOptional('xtrue',[]);
            p.parse(varargin{:});
            opts = p.Results;
            
            % Extract data from obj
            R = obj.R;
            n = size(R,1);
            a = obj.alpha;
            v = obj.v;
            
            niter = opts.maxiter;
            tol = opts.tol;
            xcur = zeros(n,1);
            xcur = xcur + v;
            
            hist = zeros(niter, 1);
            
            I = eye(n);
            
            for i = 1:niter
                A = kron(xcur, I) + kron(I, xcur);
                A = I - a/2*R*A;
                b = (1-a)*v;
                xn = A \ b;
                xn = xn / norm(xn, 1);
                
                curdiff = norm(xn - xcur, 1);
                hist(i) = curdiff;
                
                if ~isempty(opts.xtrue), hist(i) = norm(xn - opts.xtrue,inf); end
                
                if curdiff <= tol
                    break
                end
           
                xcur = xn;
            end
            
            
            hist = hist(1:i, :);
            if i == niter && curdiff > tol
                warning('did not converge');
                flag = 0;
            else
                flag = 1;
            end
            
            x = xn;
        end
        
        
        function [x, hist, flag] = newton(obj, varargin)
            % NEWTON Solve the tensorpr3 iteration using Newton's method
            
            p = inputParser;
            p.addOptional('maxiter',1e5);
            p.addOptional('tol',1e-8);
            p.addOptional('xtrue',[]);
            p.parse(varargin{:});
            opts = p.Results;
            
            % Extract data from obj
            R = obj.R;
            n = size(R,1);
            a = obj.alpha;
            v = obj.v;
            
            niter = opts.maxiter;
            tol = opts.tol;
            xcur = zeros(n,1);
            xcur = xcur + v;
            
            hist = zeros(niter, 1);
            
            I = eye(n);
            for i = 1:niter
                A = a*R*(kron(xcur, I) + kron(I, xcur)) - I;
                b = a*R*kron(xcur, xcur) - (1-a)*v; % residual
                xn = A \ b;
                xn = xn / sum(xn);
                
                curdiff = norm(xn - xcur, 1);
                hist(i) = curdiff;
                
                if ~isempty(opts.xtrue), hist(i) = norm(xn - opts.xtrue,inf); end
                
                if curdiff <= tol
                    break
                end
           
                xcur = xn;
            end
            
            hist = hist(1:i, :);
            if i == niter && curdiff > tol
                warning('did not converge');
                flag = 0;
            else
                flag = 1;
            end
            
            x = xn;
        end
    
        function [x, hist, flag] = inner_outer(obj, varargin)
            % inner_outer method
            p = inputParser;
            p.addOptional('maxiter',1e5);
            p.addOptional('tol',1e-8);
            p.addOptional('xtrue',[]);
            p.parse(varargin{:});
            opts = p.Results;
            niter = opts.maxiter;
            tol = opts.tol;
            
            % Extract data from obj
            R = obj.R;
            n = size(R,1);
            a = obj.alpha;
            v = obj.v;
            
            Rt = a*R + (1-a)*v*ones(1, n^2);
            Rt = normout(Rt);
            at = a / 2;
            xt = v;
            hist = zeros(niter, 1);
            for i = 1:niter
                Tr = tensorpr3(Rt, at, xt);
                xt2 = Tr.solve();
                curdiff = norm(xt - xt2, 1);
                hist(i) = curdiff;
                if ~isempty(opts.xtrue)
                    hist(i) = norm(xt2 - opts.xtrue,inf);
                end
                
                % check for termination
                if curdiff <= tol
                    break;
                end
                % switch solutions
                xt = xt2;                
            end
            hist = hist(1:i,:);
            if i == niter && curdiff > tol
                warning('did not converge');
                flag = 0;
            else
                flag = 1;
            end
            
            x = xt2;
        end
        
        function [P,MR] = markov(obj)
            % Return the tensors and matrices for the modified Markov chain 
            n = size(obj.R,1);
            e = ones(n,1);
            MR = obj.alpha * obj.R + (1-obj.alpha)*(obj.v * kron(e',e'));
            P = reshape(obj.R,n,n,n);
        end
    
        function P = markov2(obj)
            % Return the transition matrix for the 2nd order Markov chain. 
            n = size(obj.R, 1);
            e = ones(n,1);
            MR = obj.alpha * obj.R + (1-obj.alpha)*(obj.v * kron(e',e'));
            n = size(MR, 1);
            P = zeros(n^2, n^2);
            for i = 1:n     % group i
                tmp = zeros(n^2, n);
                for j = 1:n % column j
                    ej = zeros(n, 1);
                    ej(j) = 1;
                    tmp(:, j) = kron(ej, MR(:, (i-1)*n + j));
                end
                P(:, (i-1)*n +1: i*n) = tmp;
            end
        end
    
    end
end