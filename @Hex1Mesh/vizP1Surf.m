function vizP1Surf(T)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if isempty(T.Surface)
    error('Surface must exist!')
end

surfX = T.SurfacePoints;

xx = reshape(surfX(:,1),3,[]);
yy = reshape(surfX(:,2),3,[]);
zz = reshape(surfX(:,3),3,[]);
xfigure;
patch(xx,yy,zz,'r');
axis equal
view(130,14)

ntri = length(surfX)/3;
nCutEle = T.SurfaceInfo.NCutElements;

title([num2str(ntri),' triangles on ',num2str(nCutEle),' cut Hex1 elements'])

end

