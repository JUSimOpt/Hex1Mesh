function [fi, fix, fiy, fiz, vol] = baseHexP1(T,iel,X)
%BASEHEXP1 Basis functions of Hex1 elements
% Computes base functions and derivatives of HEX element in one or several
% points
% [fi, fix, fiy, fiz, vol] = baseHexP1(T,iel,X)
% T is the HEX mesh class
% iel is the element index number
% X is a set of coordinates at which to evaluate the basefunctions.

iv = T.Connectivity(iel,:);
% I = [1     2     4     3     5     8     6     7];
% xc = T.XC(iv(I)); yc = T.YC(iv(I)); zc = T.ZC(iv(I));
xc = T.XC(iv); yc = T.YC(iv); zc = T.ZC(iv);


A = [1, xc(1),yc(1),zc(1),xc(1)*yc(1),yc(1)*zc(1),zc(1)*xc(1),xc(1)*yc(1)*zc(1);...
    1, xc(2),yc(2),zc(2),xc(2)*yc(2),yc(2)*zc(2),zc(2)*xc(2),xc(2)*yc(2)*zc(2);...
    1, xc(3),yc(3),zc(3),xc(3)*yc(3),yc(3)*zc(3),zc(3)*xc(3),xc(3)*yc(3)*zc(3);...
    1, xc(4),yc(4),zc(4),xc(4)*yc(4),yc(4)*zc(4),zc(4)*xc(4),xc(4)*yc(4)*zc(4);...
    1, xc(5),yc(5),zc(5),xc(5)*yc(5),yc(5)*zc(5),zc(5)*xc(5),xc(5)*yc(5)*zc(5);...
    1, xc(6),yc(6),zc(6),xc(6)*yc(6),yc(6)*zc(6),zc(6)*xc(6),xc(6)*yc(6)*zc(6);...
    1, xc(7),yc(7),zc(7),xc(7)*yc(7),yc(7)*zc(7),zc(7)*xc(7),xc(7)*yc(7)*zc(7);...
    1, xc(8),yc(8),zc(8),xc(8)*yc(8),yc(8)*zc(8),zc(8)*xc(8),xc(8)*yc(8)*zc(8)];

fim = (A\eye(8));


np = size(X,1);
o1 = ones(1,np);
z1 = zeros(1,np);

x = X(:,1); y = X(:,2); z = X(:,3);
x = x(:)'; y = y(:)'; z = z(:)';
fi = fim'*[o1;x;y;z;x.*y;y.*z;z.*x;x.*y.*z];

vol = (xc(3)-xc(1))*(yc(2)-yc(1))*(zc(5)-zc(1));

fix = fim'*[z1;o1;z1;z1;y;z1;z;y.*z];
fiy = fim'*[z1;z1;o1;z1;x;z;z1;x.*z];
fiz = fim'*[z1;z1;z1;o1;z1;y;x;x.*y];

