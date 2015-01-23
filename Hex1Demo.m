clear
close all
clc



x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 4;
nye = 4;
nze = 4;

x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 3;
nye = 3;
nze = 3;

disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

hv = mesh1.vizMesh('ElementNumbers','NodeNumbers');


mesh1.mesh

nodes = mesh1.Connectivity;
xnod = mesh1.xnod; ynod = mesh1.ynod; znod = mesh1.znod;
X= [xnod,ynod,znod];
mh = mesh1.mesh

tic
N1 = mesh1.Neighbors('Structured');
toc
tic
N2 = mesh1.Neighbors('Naive');
toc

X = [mesh1.xnod,mesh1.ynod,mesh1.znod];
nodes = mesh1.Connectivity;

%%
tic
neighs = ComputeHexNeighborsNaive(X,nodes);
toc
%%


