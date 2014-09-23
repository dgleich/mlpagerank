function rval = strong_connectivity_filter(T)

Z = sum(T,1);
if all(Z(:)>0)
    M = sum(T,3);
    rval = max(components(sparse(M)))==1;
else
    rval = 0;
end
