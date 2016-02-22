clear
close all
clc


x0 = -1;
x1 = 1;
y0 = -1;
y1 = 1;
z0 = -1;
z1 = 1;
ne = 2;
nxe = ne;nye = ne;nze = ne;

disp('creating mesh...')
tic
H = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

hv = H.vizMesh('ElementNumbers','NodeNumbers');
H.RefineLocal([1,8])
hv = H.vizMesh('ElementNumbers','NodeNumbers');
return
N1 = H.Neighbors('Structured');


xnod = H.XC;
ynod = H.YC;
znod = H.ZC;
%% Surface function
% Create the surface function
R = 0.89;
xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
% xc = 0; yc = 0; zc = 0;
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;
surfaceFunction = @(x,y) ((x - xc).^2+(y - yc).^2).^.5-R;

% phi = surfaceFunction(xnod, ynod, znod);
phi = surfaceFunction(xnod, ynod);


%% CutP1 Surface
% disp('CutP1')
% tic
% H.CutP1(phi);
% toc

%% Triangulate P1
% disp('Triangulate P1')
% tic
% [tri,surfX] = H.TriangulateP1;
% toc
% h2 = H.vizP1Surf()

%% CutP2 Surface
surfh = H.CutP2(phi,0);
%% Viz P2 Surface
disp('Viz P2 surface')
tic
hv = H.vizP2Surf('FaceColor','c','EdgeColor','k','nP2Ele',5);
toc



%% Basis functions and numerical integration

















