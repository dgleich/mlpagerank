function S = combine_structs(varargin)
% COMBINE_STRUCTS Create a combination of structs with members from each
S = varargin{1};
assert(isstruct(S));
for i=2:numel(varargin)
    Si = varargin{i};
    assert(isstruct(Si));
    fn = fieldnames(Si);
    for j=1:numel(fn)
        S.(fn{j}) = Si.(fn{j});
    end
end