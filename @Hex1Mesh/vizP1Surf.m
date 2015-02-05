function hp = vizP1Surf(T)
%VIZP1SURF Visualize P1 surface
%   hp = vizP1Surf(T)

if isempty(T.SurfaceP1)
    error('Surface must exist!')
end

if ~isempty(T.SurfaceP1Triangulation)
    hp.fig = xfigure; axis equal; hold on;
    FV.Vertices = T.SurfaceP1Points;
    FV.Faces = T.SurfaceP1Triangulation;
    % shading interp
    hp.light = light;
    hp.patch = patch(FV,'FaceColor','c','FaceLighting','gouraud');
    [az,el]=view(136,16);
    hp.az = az;
    hp.el = el;

else
    surfX = T.SurfaceP1Points;
    xx = reshape(surfX(:,1),3,[]);
    yy = reshape(surfX(:,2),3,[]);
    zz = reshape(surfX(:,3),3,[]);
    hp.fig = xfigure; axis equal; hold on;
    hp.patch = patch(xx,yy,zz,'r');
    [az,el]=view(130,14);
    hp.az = az;
    hp.el = el;
end

ntri = length(T.SurfaceP1Points)/3;
nCutEle = T.SurfaceP1Info.NCutElements;

title([num2str(ntri),' triangles on ',num2str(nCutEle),' cut Hex1 elements'])

end

