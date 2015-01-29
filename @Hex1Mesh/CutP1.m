function [surfX, surfh] = CutP1(T, phi)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nodes = T.Connectivity;
xnod = T.XC;
ynod = T.YC;
znod = T.ZC;

SurfEle = find(sum(phi(nodes)>0,2) > 0 & sum(phi(nodes)>0,2)< 8);

nprec = 1e-8;
isequalAbs = @(x,y) ( abs(x-y) <= nprec );
isodd = @(x) mod(x,2);
isequalPoint = @(a,b) (sum(isequalAbs(a, b))==3);

nsurfEle = length(SurfEle);

% preallocate space for coords
% Maximum size is two triangles per HEX
% Set all to nan rather than zero to be able to find empty elements later
surfX = ones(nsurfEle*2,3)*nan;
surfh(nsurfEle*2).iel = [];
nTriEle = 1;% number of triangles found

for iel = SurfEle'
    %% Compute cut points
    iv = nodes(iel,:);
    faces = T.Element(iel).faces;
    Pe = NaN(2*6,3);
    ed = 1;
    normal = zeros(12,3);
    for j = 1:size(faces,1)
        ivv = faces(j,:);
        xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
        XC =[xc,yc,zc];
        q  = phi(faces(j,:));
        %         for in = 1:4
        %             text(xc(in),yc(in),zc(in), [num2str(ivv(in)),', ',num2str(q(in))] )
        %         end
        %         axis equal tight;
        
        locEdges = [ivv(1),ivv(2);...
            ivv(2),ivv(3);...
            ivv(3),ivv(4);...
            ivv(4),ivv(1)];
        %         text(mean(xc),mean(yc),mean(zc),num2str(j))
        
        for ied = 1:4
            if ied == 4
                ui = [q(4),q(1)];
                xi = [XC(4,:);XC(1,:)];
            else
                ui = [q(ied),q(ied+1)];
                xi = [XC(ied,:);XC(ied+1,:)];
            end
            if sign(ui(1)) ~= sign(ui(2))
                %                disp('edge cut')
                %                locEdges(ied,:)
                [c,k] = sort(ui);
                Xi = xi(k,:);
                ip = Xi(1,:)+(Xi(2,:)-Xi(1,:)).*((0-c(1))/(c(2)-c(1)));
                Pe(ed,:) = ip;
                n = Xi(2,:)-Xi(1,:);
                normal(ed,:) = n;
                %                 plot3(ip(1),ip(2),ip(3),'b*')
                ed = ed+1;
            end
            
        end
        
    end
    
    %% Convert to triangle elements
    normal = sum(normal,1);
    normal = normal/norm(normal);
    
    Xe = unique(Pe,'rows','stable');
    Xe = Xe(~any(isnan(Xe),2),:);
    
    
    if size(Xe,1) == 3
        %% Triangle
        t1 = [1,2,3];
        
        Xp = Xe;
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t1 = [1,3,2];
            n = -n;
        end
        
        %         patch(Xe(t1,1),Xe(t1,2),Xe(t1,3),'c','EdgeColor','k')
        %         quiver3(mean(Xe(t1,1)),mean(Xe(t1,2)),mean(Xe(t1,3)),n(1),n(2),n(3),0.1,'Color','y')
        
        % Add to surfh
        surfh(nTriEle).iel  = iel;
        surfh(nTriEle).Xe    = Xe(t1,:);
        surfh(nTriEle).xp = Xe(t1,1);
        surfh(nTriEle).yp = Xe(t1,2);
        surfh(nTriEle).zp = Xe(t1,3);
        %         surfh(nTriEle).nei = neighs(iel,:);
        surfh(nTriEle).faceNormal = n;
        
        %Assembly
        lo = nTriEle*3-2;
        up = nTriEle*3;
        surfX(lo:up,1:3) = Xe(t1,:);
        nTriEle = nTriEle +1;
        
    elseif size(Xe,1) == 4
        %% Quadliteral
        % Viz Polygon normal
        %         for i = 1:4
        %             text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
        %         end
        %         quiver3(mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3)),normal(1),normal(2),normal(3),0.1,'Color','k')
        
        % Rotate the polygon down to the x-y plane
        origin = Xe(1,:);
        localz = cross(Xe(2,:)-origin, Xe(3,:)-origin);
        unitz = localz/norm(localz,2);
        %calculate local x vector in plane
        localx = Xe(2,:)-origin;
        unitx = localx/norm(localx,2);
        %calculate local y
        localy = cross(localz, localx);
        unity = localy/norm(localy,2);
        TM = [localx(:), localy(:), localz(:), origin(:); 0 0 0 1];
        CM = [Xe, ones(4,1)];
        Xe2D = TM \ CM';
        Xe2D = Xe2D(1:3,:)';
        % Measure the curvature
        sx = max(Xe2D(:,1))-min(Xe2D(:,1));
        sy = max(Xe2D(:,2))-min(Xe2D(:,2));
        sz = max(Xe2D(:,3))-min(Xe2D(:,3));
        curvature = (sz / ((sx*sy)^(1/2)))*100;
        if curvature > 1
            warning('excessive curvature in cut elements! Increase number of elements')
        end
        % New order of polygon points
        k = convhull(Xe2D(:,1),Xe2D(:,2));
        k1 = k(1:end-1);
        Xe = Xe(k1,:);
        
        % Polygon edge color
        %         plot3(Xe([1:4,1],1),Xe([1:4,1],2),Xe([1:4,1],3),'-k')
        
        % Split into triangles
        t1 = [1,3,2];
        t2 = [1,4,3];
        
        % triangle 1
        Xp = Xe(t1,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t1 = [1,2,3];
            n = -n;
        end
        n1 = n;
        %         patch(Xe(t1,1),Xe(t1,2),Xe(t1,3),'b','EdgeColor','none')
        %         quiver3(mean(Xe(t1,1)),mean(Xe(t1,2)),mean(Xe(t1,3)),n(1),n(2),n(3),0.1,'Color','g')
        
        % triangle 2
        Xp = Xe(t2,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t2 = [1,3,4];
            n = -n;
        end
        n2 = n;
        %         patch(Xe(t2,1),Xe(t2,2),Xe(t2,3),'b','EdgeColor','none')
        %         quiver3(mean(Xe(t2,1)),mean(Xe(t2,2)),mean(Xe(t2,3)),n(1),n(2),n(3),0.1,'Color','g')
        
        % Add to surfh
        surfh(nTriEle).iel  = iel;
        surfh(nTriEle).Xe    = Xe(t1,:);
        surfh(nTriEle).xp = Xe(t1,1);
        surfh(nTriEle).yp = Xe(t1,2);
        surfh(nTriEle).zp = Xe(t1,3);
        %         surfh(nTriEle).nei = neighs(iel,:);
        surfh(nTriEle).faceNormal = n1;
        
        surfh(nTriEle+1).iel  = iel;
        surfh(nTriEle+1).Xe    = Xe(t1,:);
        surfh(nTriEle+1).xp = Xe(t1,1);
        surfh(nTriEle+1).yp = Xe(t1,2);
        surfh(nTriEle+1).zp = Xe(t1,3);
        %         surfh(nTriEle+1).nei = neighs(iel,:);
        surfh(nTriEle+1).faceNormal = n1;
        
        % Assembly
        lo = nTriEle*3-2;
        up = nTriEle*3;
        surfX(lo:up,1:3) = Xe(t1,:);
        
        lo = (nTriEle+1)*3-2;
        up = (nTriEle+1)*3;
        surfX(lo:up,1:3) = Xe(t2,:);
        
        nTriEle = nTriEle +2;
        
    elseif size(Xe,1) == 5
        %% Pentagon
        % Rotate the polygon down to the x-y plane
        origin = Xe(1,:);
        localz = cross(Xe(2,:)-origin, Xe(3,:)-origin);
        unitz = localz/norm(localz,2);
        %calculate local x vector in plane
        localx = Xe(2,:)-origin;
        unitx = localx/norm(localx,2);
        %calculate local y
        localy = cross(localz, localx);
        unity = localy/norm(localy,2);
        TM = [localx(:), localy(:), localz(:), origin(:); 0 0 0 1];
        CM = [Xe, ones(5,1)];
        Xe2D = TM \ CM';
        Xe2D = Xe2D(1:3,:)';
        % Measure the curvature
        sx = max(Xe2D(:,1))-min(Xe2D(:,1));
        sy = max(Xe2D(:,2))-min(Xe2D(:,2));
        sz = max(Xe2D(:,3))-min(Xe2D(:,3));
        curvature = (sz / ((sx*sy)^(1/2)))*100;
        if curvature > 1
            warning('excessive curvature in cut elements! Increase number of elements')
        end
        
        % New order of polygon points
        k = convhull(Xe2D(:,1),Xe2D(:,2));
        k1 = k(1:end-1);
        Xe = Xe(k1,:);
        
        % Viz polygon edge
        %         plot3(Xe([1:5,1],1),Xe([1:5,1],2),Xe([1:5,1],3),'-k')
        %         quiver3(mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3)),normal(1),normal(2),normal(3),0.1,'Color','k')
        
        % Split into triangles
        t1 = [1,2,3];
        t2 = [1,3,4];
        t3 = [1,4,5];
        
        % triangle 1
        Xp = Xe(t1,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t1 = [1,3,2];
            n = -n;
        end
        n1 = n;
        %         patch(Xe(t1,1),Xe(t1,2),Xe(t1,3),'r','EdgeColor','none')
        %         quiver3(mean(Xe(t1,1)),mean(Xe(t1,2)),mean(Xe(t1,3)),n(1),n(2),n(3),0.1,'Color','b')
        
        % triangle 2
        Xp = Xe(t2,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t2 = [1,4,3];
            n = -n;
        end
        n2 = n;
        %         patch(Xe(t2,1),Xe(t2,2),Xe(t2,3),'r','EdgeColor','none')
        %         quiver3(mean(Xe(t2,1)),mean(Xe(t2,2)),mean(Xe(t2,3)),n(1),n(2),n(3),0.1,'Color','b')
        
        % triangle 3
        Xp = Xe(t3,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t3 = [1,5,4];
            n = -n;
        end
        n3 = n;
        %         patch(Xe(t3,1),Xe(t3,2),Xe(t3,3),'r','EdgeColor','none')
        %         quiver3(mean(Xe(t3,1)),mean(Xe(t3,2)),mean(Xe(t3,3)),n(1),n(2),n(3),0.1,'Color','b')
        
        % Add to surfh
        surfh(nTriEle).iel = iel;
        surfh(nTriEle).Xe = Xe(t1,:);
        surfh(nTriEle).xp = Xe(t1,1);
        surfh(nTriEle).yp = Xe(t1,2);
        surfh(nTriEle).zp = Xe(t1,3);
        %         surfh(nTriEle).nei = neighs(iel,:);
        surfh(nTriEle).faceNormal = n1;
        
        surfh(nTriEle+1).iel = iel;
        surfh(nTriEle+1).Xe = Xe(t2,:);
        surfh(nTriEle+1).xp = Xe(t2,1);
        surfh(nTriEle+1).yp = Xe(t2,2);
        surfh(nTriEle+1).zp = Xe(t2,3);
        %         surfh(nTriEle+1).nei = neighs(iel,:);
        surfh(nTriEle+1).faceNormal = n2;
        
        surfh(nTriEle+2).iel = iel;
        surfh(nTriEle+2).Xe = Xe(t3,:);
        surfh(nTriEle+2).xp = Xe(t3,1);
        surfh(nTriEle+2).yp = Xe(t3,2);
        surfh(nTriEle+2).zp = Xe(t3,3);
        %         surfh(nTriEle+2).nei = neighs(iel,:);
        surfh(nTriEle+2).faceNormal = n3;
        
        % Assembly
        lo = nTriEle*3-2;
        up = nTriEle*3;
        surfX(lo:up,1:3) = Xe(t1,:);
        
        lo = (nTriEle+1)*3-2;
        up = (nTriEle+1)*3;
        surfX(lo:up,1:3) = Xe(t2,:);
        
        lo = (nTriEle+2)*3-2;
        up = (nTriEle+2)*3;
        surfX(lo:up,1:3) = Xe(t3,:);
        
        nTriEle = nTriEle +3;
        
    elseif size(Xe,1) == 6
        %% Hexagon
        % Rotate the polygon down to the x-y plane
        origin = Xe(1,:);
        localz = cross(Xe(2,:)-origin, Xe(3,:)-origin);
        unitz = localz/norm(localz,2);
        %calculate local x vector in plane
        localx = Xe(2,:)-origin;
        unitx = localx/norm(localx,2);
        %calculate local y
        localy = cross(localz, localx);
        unity = localy/norm(localy,2);
        TM = [localx(:), localy(:), localz(:), origin(:); 0 0 0 1];
        CM = [Xe, ones(6,1)];
        Xe2D = TM \ CM';
        Xe2D = Xe2D(1:3,:)';
        % Measure the curvature
        sx = max(Xe2D(:,1))-min(Xe2D(:,1));
        sy = max(Xe2D(:,2))-min(Xe2D(:,2));
        sz = max(Xe2D(:,3))-min(Xe2D(:,3));
        curvature = (sz / ((sx*sy)^(1/2)))*100;
        if curvature > 1
            warning('excessive curvature in cut elements! Increase number of elements')
        end
        
        % New order of polygon points
        k = convhull(Xe2D(:,1),Xe2D(:,2));
        k1 = k(1:end-1);
        Xe = Xe(k1,:);
        
        % Viz polygon edge
        %         plot3(Xe([1:6,1],1),Xe([1:6,1],2),Xe([1:6,1],3),'-k')
        
        %         patch(Xe(:,1),Xe(:,2),Xe(:,3),'g')
        %         quiver3(mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3)),normal(1),normal(2),normal(3),0.1,'Color','k')
        %         for i = 1:6
        %            text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
        %         end
        %         pause
        
        % Split into triangles
        t1 = [1,2,3];
        t2 = [1,3,4];
        t3 = [1,4,5];
        t4 = [1,5,6];
        
        % triangle 1
        Xp = Xe(t1,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t1 = [1,3,2];
            n = -n;
        end
        n1 = n;
        %         patch(Xe(t1,1),Xe(t1,2),Xe(t1,3),'m','EdgeColor','none')
        %         quiver3(mean(Xe(t1,1)),mean(Xe(t1,2)),mean(Xe(t1,3)),n(1),n(2),n(3),0.1,'Color','k')
        
        % triangle 2
        Xp = Xe(t2,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t2 = [1,4,3];
            n = -n;
        end
        n2 = n;
        %         patch(Xe(t2,1),Xe(t2,2),Xe(t2,3),'m','EdgeColor','none')
        %         quiver3(mean(Xe(t2,1)),mean(Xe(t2,2)),mean(Xe(t2,3)),n(1),n(2),n(3),0.1,'Color','k')
        
        % triangle 3
        Xp = Xe(t3,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t3 = [1,5,4];
            n = -n;
        end
        n3 = n;
        %         patch(Xe(t3,1),Xe(t3,2),Xe(t3,3),'m','EdgeColor','none')
        %         quiver3(mean(Xe(t3,1)),mean(Xe(t3,2)),mean(Xe(t3,3)),n(1),n(2),n(3),0.1,'Color','k')
        
        % triangle 4
        Xp = Xe(t4,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/ (n(1)^2+n(2)^2+n(3)^2)^0.5;
        if sign(n) ~= sign(normal)
            t4 = [1,6,5];
            n = -n;
        end
        n4 = n;
        %         patch(Xe(t4,1),Xe(t4,2),Xe(t4,3),'m','EdgeColor','none')
        %         quiver3(mean(Xe(t4,1)),mean(Xe(t4,2)),mean(Xe(t4,3)),n(1),n(2),n(3),0.1,'Color','k')
        
        % Add to surfh
        surfh(nTriEle).iel = iel;
        surfh(nTriEle).Xe = Xe(t1,:);
        surfh(nTriEle).xp = Xe(t1,1);
        surfh(nTriEle).yp = Xe(t1,2);
        surfh(nTriEle).zp = Xe(t1,3);
        %         surfh(nTriEle).nei = neighs(iel,:);
        surfh(nTriEle).faceNormal = n1;
        
        surfh(nTriEle+1).iel = iel;
        surfh(nTriEle+1).Xe = Xe(t2,:);
        surfh(nTriEle+1).xp = Xe(t2,1);
        surfh(nTriEle+1).yp = Xe(t2,2);
        surfh(nTriEle+1).zp = Xe(t2,3);
        %         surfh(nTriEle+1).nei = neighs(iel,:);
        surfh(nTriEle+1).faceNormal = n2;
        
        surfh(nTriEle+2).iel = iel;
        surfh(nTriEle+2).Xe = Xe(t3,:);
        surfh(nTriEle+2).xp = Xe(t3,1);
        surfh(nTriEle+2).yp = Xe(t3,2);
        surfh(nTriEle+2).zp = Xe(t3,3);
        %         surfh(nTriEle+2).nei = neighs(iel,:);
        surfh(nTriEle+2).faceNormal = n3;
        
        surfh(nTriEle+3).iel = iel;
        surfh(nTriEle+3).Xe = Xe(t4,:);
        surfh(nTriEle+3).xp = Xe(t4,1);
        surfh(nTriEle+3).yp = Xe(t4,2);
        surfh(nTriEle+3).zp = Xe(t4,3);
        %         surfh(nTriEle+3).nei = neighs(iel,:);
        surfh(nTriEle+3).faceNormal = n4;
        
        % Assembly
        lo = nTriEle*3-2;
        up = nTriEle*3;
        surfX(lo:up,1:3) = Xe(t1,:);
        
        lo = (nTriEle+1)*3-2;
        up = (nTriEle+1)*3;
        surfX(lo:up,1:3) = Xe(t2,:);
        
        lo = (nTriEle+2)*3-2;
        up = (nTriEle+2)*3;
        surfX(lo:up,1:3) = Xe(t3,:);
        
        lo = (nTriEle+3)*3-2;
        up = (nTriEle+3)*3;
        surfX(lo:up,1:3) = Xe(t4,:);
        
        nTriEle = nTriEle +4;
        
    else
        %         patch(Xe(:,1),Xe(:,2),Xe(:,3),'k')
        pause
    end
    %     drawnow
end

T.Surface = surfh;
T.SurfacePoints = surfX;
T.SurfaceInfo.CutElements = SurfEle;
T.SurfaceInfo.NCutElements = length(SurfEle);



end

