clear
close all
clc


x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
ne = 10;
nxe = ne;nye = ne;nze = ne;

disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

% hv = mesh1.vizMesh('ElementNumbers','NodeNumbers');

N1 = mesh1.Neighbors('Structured');


xnod = mesh1.XC;
ynod = mesh1.YC;
znod = mesh1.ZC;
%% Surface function
% Create the surface function
R = 0.89;
% surfaceFunction = @(x,y) ((x - x0).^2+(y - y0).^2).^.5-R;
% xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
xc = 0; yc = 0; zc = 0;
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;

phi = surfaceFunction(xnod, ynod, znod);


%% Surface
tic
[surfX, surfh] = mesh1.CutP1(phi);
CutP1Time = toc

mesh1.vizP1Surf()

%% Triangulation
