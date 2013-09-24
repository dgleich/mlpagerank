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
                v = 1/n;
            end
            
            obj.R = R;
            obj.alpha = alpha;
            obj.v = v;
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
            
            Gamma = 1 / (1+gamma);
            xcur = zeros(n,1); % this allows us to keep v = 1/n :-)
            xcur = xcur + v; 
            
            hist = zeros(niter, 1);
            
            for i=1:niter
                % TODO make this iteration better
                y = a*(R*kron(xcur, xcur)); % make sure Matlab does it 
                z = y * Gamma + Gamma*(1-sum(y))*v;
                xn = z + (1-sum(z))*xcur;
                
                curdiff = norm(xn - xcur, inf);
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
            
            x = xn;
        end
    end
end
