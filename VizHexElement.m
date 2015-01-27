function  VizHexElement(T, ele)

nodes = T.Connectivity;
xnod = T.xnod;
ynod = T.ynod;
znod = T.znod;
nele = size(nodes,1);

if strcmpi(ele,'all')
    ele = 1:nele;
end

XC = zeros(4,length(ele)*6);YC = XC;ZC = XC;
c = 1;
for iel = ele
   iv = nodes(iel,:);
   faces = T.Element(iel).faces;
   for j = 1:size(faces,1)
       ivv = faces(j,:);
       xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
       
       XC(1:4,c)=xc;
       YC(1:4,c)=yc;
       ZC(1:4,c)=zc;
       c=c+1;
   end
   
end

xfigure;
patch(XC,YC,ZC,'w','FaceColor','none');  hold on;view(23,37),axis equal
end

