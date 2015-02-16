clear
close all
clc


x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
ne = 6;
nxe = ne;nye = ne;nze = ne;

disp('creating mesh...')
tic
H = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

% hv = H.vizMesh('ElementNumbers','NodeNumbers');

N1 = H.Neighbors('Structured');


xnod = H.XC;
ynod = H.YC;
znod = H.ZC;
%% Surface function
% Create the surface function
R = 0.89;
% surfaceFunction = @(x,y) ((x - x0).^2+(y - y0).^2).^.5-R;
% xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
xc = 0; yc = 0; zc = 0;
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;

phi = surfaceFunction(xnod, ynod, znod);


%% CutP1 Surface
disp('CutP1')
tic
H.CutP1(phi);
toc

%% Triangulate P1
disp('Triangulate P1')
tic
[tri,surfX] = H.TriangulateP1;
toc
h2 = H.vizP1Surf()

%% CutP2 Surface
surfh = H.CutP2(phi);
%% Viz P2 Surface
disp('Viz P2 surface')
tic
hv = H.vizP2Surf('FaceColor','c','EdgeColor','k','nP2Ele',5);
toc



%% Basis functions and numerical integration

















