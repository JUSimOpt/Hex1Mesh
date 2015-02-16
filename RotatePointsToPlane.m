function [k1] = RotatePointsToPlane(Xe)
    origin = Xe(1,:);
    localz = cross(Xe(2,:)-origin, Xe(3,:)-origin);
    unitz = localz/norm(localz,2);
    %calculate local x vector in plane
    localx = Xe(2,:)-origin;
    unitx = localx/norm(localx,2);
    %calculate local y
    localy = cross(localz, localx);
    unity = localy/norm(localy,2);
    TM = [unitx(:), unity(:), unitz(:), origin(:); 0 0 0 1];
    CM = [Xe, ones(size(Xe,1),1)];
    Xe2D = TM \ CM';
    Xe2D = Xe2D(1:3,:)';
    % Measure the curvature
    sx = max(Xe2D(:,1))-min(Xe2D(:,1));
    sy = max(Xe2D(:,2))-min(Xe2D(:,2));
    sz = max(Xe2D(:,3))-min(Xe2D(:,3));
    curvature = (sz / mean([sx,sy]))*100;
    if curvature > 10
        warning(['excessive curvature (',num2str(curvature),') in cut elements! Increase number of elements!'])
% 
%         Xe2D
%         xfigure(55);clf; hold on;axis equal;
%         sx
%         sy
%         sz
%         sz/sx
%         putvar(Xe2D,sx,sy,sz)
%         curvature
%         plot3(Xe2D(:,1),Xe2D(:,2),Xe2D(:,3),'*k')
%         pause
    end
    
%     xfigure(55);clf; hold on;axis equal;
%     plot3(Xe2D(:,1),Xe2D(:,2),Xe2D(:,3),'*k')
    
    % New order of polygon points
    k = convhull(Xe2D(:,1),Xe2D(:,2));
    k1 = k(1:end-1);
end