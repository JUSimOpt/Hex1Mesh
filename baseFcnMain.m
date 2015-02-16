clear
clc
close all

%%
x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 1;
nye = 1;
nze = 1;

iel = 1;

disp('creating mesh...')
tic
H = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc
hv = H.vizMesh(iel,'ElementNumbers','NodeNumbers');
%     8-----7
%    /|    /|
%   6-----5 |
%   | 4...|.3
%   |/    |/ 
%   2-----1
%%
% X = rand(1,3)
X = [0,0,0;...
     0,1,0;...
     1,0,0;...
     1,1,0;...
     0,0,1;...
     0,1,1;...
     1,0,1;...
     1,1,1];

 
tic
[fi, fix, fiy, fiz, vol] = baseHexP1(H,iel,X);
toc

%%
nele = 5000;
disp(['Calling baseHex for 8 points over ',num2str(nele),' elements'])
tic
for i = 1:nele
    X = rand(8,3);
    [fi, fix, fiy, fiz, vol] = baseHexP1(H,iel,X);
end
toc
% 
% disp(['Calling baseHexCubic for 8 points over ',num2str(nele),' elements'])
% tic
% for i = 1:nele
%     X = rand(8,3);
%     [fi] = baseHexCubic(mesh1,iel,X(:,1),X(:,2),X(:,3));
% end
% toc

%%
%syms x y z
%[fi, fix, fiy, fiz, vol] = baseHex(mesh1,iel,x,y,z)




