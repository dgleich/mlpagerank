function x = subset_enumerate(cur,bins)
% SUBSET_ENUMERATOR Enumerate all subsets with a generalized counter
%
% Implement a generalized counter: set bins(i) to be the total choices for
% the ith position of the counter, then 
% next = subset_enumerator(cur,bins)
% returns another valid counter where 1 <= next(i) <= bins(i) and 
% successive enumerate will identify ALL possible values.  If bins(i)=2 for
% all i, then this just implements a length(bins)-bit binary counter, but
% it can be useful for exact subset enumeration.
%
% Example:
%   bins = 2*ones(5,1); % look at 
%   cur = subset_enumerate(bins)
%   while (cur ~= 0)
%     fprintf('[ '); fprintf('%i ', cur); fprintf(']\n');
%     cur = subset_enumerate(cur,bins);
%   end

% Copyright, David F. Gleich 2010

% History
% :2010-02-16: Initial writing

if nargin==1 
    % initialization
    bins = cur;
    x = ones(1,length(bins));
else
    % nargin == 2...
    x = cur;
    n = length(bins);
    j = n;
    while j>=1 && x(j)==bins(j)
        j=j-1; 
    end
    if j<1, x = 0; return; end % this indicates we are done!
    x(j)=x(j)+1;
    for k=(j+1):n, x(k)=1; end
end
   
