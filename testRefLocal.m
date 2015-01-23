clear
close all
clc


x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 2;
nye = 2;
nze = 2;

disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

hv = mesh1.vizMesh('ElementNumbers','NodeNumbers');


locnodes = mesh1.Connectivity;
xnod = mesh1.xnod; ynod = mesh1.ynod; znod = mesh1.znod;
X= [xnod,ynod,znod];
% mh = mesh1.mesh


N1 = mesh1.Neighbors('Structured');

%%

mesh1.RefineLocal([2,3,4])

mesh1.vizMesh('ElementNumbers','NodeNumbers');











