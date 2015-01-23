function [fi, fix, fiy, fiz, vol] = baseHex(hm,iel,x,y,z)
% computes base functions and derivatives of HEX element in one or several
% points
% hm is the HEX mesh class
% iel is the index
% x, y and z are vector points. Basefunctions are evaluated in these
% points.
iv = hm.Connectivity(iel,:);
I = [1     2     4     3     5     8     6     7];
xc = hm.xnod(iv(I)); yc = hm.ynod(iv(I)); zc = hm.znod(iv(I));
% xc = hm.xnod(iv); yc = hm.ynod(iv); zc = hm.znod(iv);

A = [1, xc(1),yc(1),zc(1),xc(1)*yc(1),yc(1)*zc(1),zc(1)*xc(1),xc(1)*yc(1)*zc(1);...
    1, xc(2),yc(2),zc(2),xc(2)*yc(2),yc(2)*zc(2),zc(2)*xc(2),xc(2)*yc(2)*zc(2);...
    1, xc(3),yc(3),zc(3),xc(3)*yc(3),yc(3)*zc(3),zc(3)*xc(3),xc(3)*yc(3)*zc(3);...
    1, xc(4),yc(4),zc(4),xc(4)*yc(4),yc(4)*zc(4),zc(4)*xc(4),xc(4)*yc(4)*zc(4);...
    1, xc(5),yc(5),zc(5),xc(5)*yc(5),yc(5)*zc(5),zc(5)*xc(5),xc(5)*yc(5)*zc(5);...
    1, xc(6),yc(6),zc(6),xc(6)*yc(6),yc(6)*zc(6),zc(6)*xc(6),xc(6)*yc(6)*zc(6);...
    1, xc(7),yc(7),zc(7),xc(7)*yc(7),yc(7)*zc(7),zc(7)*xc(7),xc(7)*yc(7)*zc(7);...
    1, xc(8),yc(8),zc(8),xc(8)*yc(8),yc(8)*zc(8),zc(8)*xc(8),xc(8)*yc(8)*zc(8)];

fim = (A\eye(8));


nx = length(x);
o1 = ones(1,nx);
z1 = zeros(1,nx);

x = x(:)'; y = y(:)'; z = z(:)';
fi = fim'*[o1;x;y;z;x.*y;y.*z;z.*x;x.*y.*z];

vol = (xc(3)-xc(1))*(yc(2)-yc(1))*(zc(5)-zc(1));

fix = fim'*[z1;o1;z1;z1;y;z1;z1;y.*z];
fiy = fim'*[z1;z1;o1;z1;x;z;z1;x.*z];
fiz = fim'*[z1;z1;z1;o1;z1;y;x;x.*y];

