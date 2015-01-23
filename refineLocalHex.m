function [] = refineLocalHex(mesh,nele)
%REFINELOCALHEX Summary of this function goes here
%   Detailed explanation goes here
T = mesh;
xnod = T.xnod;
ynod = T.ynod;
znod = T.znod;

nele = [1,3] %What element numbers to refine

% Number of nodes, keeps track of the latest node number
nmax = length(xnod);


end

