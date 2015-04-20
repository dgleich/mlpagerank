function search_tensors(n,f,name,varargin)
% SEARCH_TENSORS Exhaustively search over binary hypermatrices for a property
%
% search_tensors(n,f,name) starts a search over tensors of size n-by-n-by-n
% tensors, or hypermatrices, looking for one where f(T) is 1. This search
% is restartable, so closing Matlab isn't a problem. Also, whenever we
% find one, we output it to a file with name 'name-found.mat'.
% Periodically, we also output a restart file with name 'name-restart.mat'
% And we always update the restart file after we find something. This is
% designed to search for tensors with some type of rare property.
%
% search_tensors([n1,...,nm],...) searches over n1-by-n2-by-...-by-nm
% tensors or hypermatrices.
%
% Options:
%    'checkpoint' : save a result every k steps (default 100000)
%    'filter' : provide a function to quickly filter out a problem
%    'max' : the largest integer considered, the default is 1, so we look
%       for binary tensors, if max=2, then we look at all tensors with 
%       non-negative integers 0, 1, 2.
%
% Simple example:
%   f = @(T) sum(T(:))>=7
%   search_tensors(2,f,'dense');
%
% Example with 4-dimensional tensor
%    search_tensors([2,2,2,2],@(T) sum(T(:))>=15,'dense');
%
% Complex example:
%   f = @(T) ~get_outputs(@() tensorpr3(normout(reshape(T,2,4)')',...
%               0.99).power,3); % get the flag output
%   search_tensors(2,f,'two-did-not-converge', ...
%       'filter',@(T) all(sum(reshape(T,2,4)) > 0));

checkname = [name '-check.mat'];
foundname = [name '-found.mat'];

p = inputParser;
p.addOptional('randomize',false);
p.addOptional('checkpoint',100000);
p.addOptional('filter', @(T) true);
p.addOptional('max',1);
p.addOptional('useperm',false);
p.addOptional('reverse',false);
p.addOptional('startnz',0);
p.parse(varargin{:});
opts = p.Results;

if numel(n)==1
    shape = {n,n,n};
else
    shape = num2cell(n);
end

if exist(checkname,'file')
    load(checkname)
else
    N = prod(cell2mat(shape));
    bins = (opts.max+1)*ones(N,1);
    cur = subset_enumerate(bins);
    iter = 1;
    found = [];
    nfound = 0;
end

if opts.useperm
    perm=randperm(N);
    %iperm(perm)=1:N;
else
    perm = 1:N;
    %iperm = 1:N;
end

reverse = opts.reverse;

if ~isempty(opts.startnz)
    cur(end-opts.startnz+1:end) = 2;
end
    

while cur ~= 0
    force_check = 0;
    T = reshape(cur(perm),shape{:}) - 1;
    if reverse
        T = opts.max - T;
    end
    if opts.filter(T)
        if f(T)
            force_check = 1;
            nfound = nfound + 1;
            S = [];
            S.(sprintf('T%i',nfound)) = T;
            if exist(foundname,'file')
                save(foundname,'-struct','S','-append');
            else
                save(foundname,'-struct','S');
            end
            fprintf('%12i   found tensor %i that passed\n',iter, nfound);
        end
    end
    iter = iter + 1;
    cur = subset_enumerate(cur, bins);
    
    if mod(iter,opts.checkpoint) == 0 || force_check
        save(checkname,'found','iter','cur','bins','nfound','perm','reverse');
    end
end
delete(checkname);