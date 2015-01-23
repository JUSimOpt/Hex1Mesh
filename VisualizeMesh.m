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


% N2 = mesh1.Neighbors('Naive');


X = [mesh1.xnod,mesh1.ynod,mesh1.znod];
nodes = mesh1.Connectivity;

T = mesh1;

iv = T.Connectivity(1,:); % local node numbers
xc = T.xnod(iv); % local coordinates
yc = T.ynod(iv);
zc = T.znod(iv);

xm = mean(xc);ym = mean(yc);zm = mean(zc); %mid points of element

edges = T.Element(1).edges; % Edge list of current element.
faces = T.Element(1).faces; % Face list of current element.
eind = [1,4,2,3,5,8,6,7,9,10,12,11];
for i = 1:12
    ei = edges(eind(i),:);
    xpi = mean(X(ei,1));
    ypi = mean(X(ei,2));
    zpi = mean(X(ei,3));
    text(xpi,ypi,zpi,[num2str(eind(i))],'BackgroundColor','Cyan')
end
faceind = [1,3,4,6,5,2]
for i = 1:6
    
    xpi = mean(X(faces(faceind(i),:),1));
    ypi = mean(X(faces(faceind(i),:),2));
    zpi = mean(X(faces(faceind(i),:),3));

    text(xpi,ypi,zpi,[num2str(faceind(i))],'BackgroundColor','r')
end
XNeind = [1,2,4,5,6,8,12,14,15,16,18,19];
XNMind = 10;
XNfind = [3,7,9,11,13,17];
