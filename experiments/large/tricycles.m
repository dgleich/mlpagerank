function [T] = tricycles(A)
% TRICYCLES Return the directed three-cycles in a graph
% T = tricycles(A) returns an array 
%   T = [ti,tj,tk]  where ti(1),tj(1),tk(1) gives one directed 3-cycle.
% 

[i,j] = find(A);
n = size(A,1);

T = tricycles_mex(n,uint32(i),uint32(j)); % take in edges