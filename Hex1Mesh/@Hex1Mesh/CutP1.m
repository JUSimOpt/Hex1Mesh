function [surfh,T] = CutP1(T, phi, level)

    nodes = T.Connectivity;
    xnod = T.XC;
    ynod = T.YC;
    znod = T.ZC;
    
    
    SurfEle = find(sum(phi(nodes)>0,2) > 0 & sum(phi(nodes)>0,2)< 8);


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
%         % Viz element
%         xfigure(1); hold on; axis equal;
%         ele = iel;
%         fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
%         patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
        
        edges = T.Element(iel).edges;
        %% Viz edges
%         ed = unique(edges);
%         for j = 1:8
%             text(xnod(ed(j)),ynod(ed(j)),znod(ed(j)),num2str(ed(j)),'BackgroundColor','w')
%         end
        

        %% Loop over all edges
        P = NaN(12,3);
        normal = P;
        for j = 1:size(edges,1)
            ivv = edges(j,:);
            xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
            XC =[xc,yc,zc];
            q  = phi(edges(j,:));
            if sign(q(1)) ~= sign(q(2))
                [u,k] = sort(q);
                Xi = XC(k,:);
                n = Xi(2,:)-Xi(1,:);
                normal(j,:) = n/norm(n);
                P(j,:) = Xi(1,:)+(Xi(2,:)-Xi(1,:)).*((level-u(1))/(u(2)-u(1)));
            end
        end
        
%         indP = find(~any(isnan(P),2));
        P = P(~any(isnan(P),2),:);
        
        normal = mean(normal(all(~isnan(normal),2),:),1);
        normal = normal/norm(normal);
        
%         quiver3(mean(P(:,1)),mean(P(:,2)),mean(P(:,3)),normal(1),normal(2),normal(3),0.1,'Color','b')
%         for i = 1:8
%            text(xnod(iv(i)),ynod(iv(i)),znod(iv(i)),num2str(phi(iv(i))),'BackgroundColor','w') 
%         end
%         plot3(P(:,1),P(:,2),P(:,3),'b*')
%         
% pause
        
        
        %% Extract triangles
        if size(P,1) == 3
            %% Triangle
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,1);
        
        elseif size(P,1) == 4
            %% Quadliteral
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,2);
                
            

        elseif size(P,1) == 5
            %% Pentagon
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,3);
                
            

        elseif size(P,1) == 6
            %% Hexagon
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,4);

        else
            patch(P(:,1),P(:,2),P(:,3),'k')
            error('Unknown element')
        end
        

    end
    
    
    nele = length([surfh.iel]);
    surfh(nele+1:end) = [];
    
    T.SurfaceP1 = surfh;
    T.SurfaceP1Points = surfX;
    T.SurfaceP1Info.CutElements = SurfEle;
    T.SurfaceP1Info.NCutElements = length(SurfEle);
    
end
    
    function [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,ntri)
    % Rotate the polygon down to the x-y plane
    k1= RotatePointsToPlane(P);

    % New order of polygon points
    Xe = P(k1,:);
    np = size(Xe,1);

%     for i = 1:4
%         text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
%     end
%     plot3(Xe([1:np,1],1),Xe([1:np,1],2),Xe([1:np,1],3),'-k')

    % Split into triangles
    tt = [ones(ntri,1),(2:ntri+1)',(3:ntri+2)'];

    for itri = 1:ntri
        % triangle itri
        ti = tt(itri,:);
        Xp = Xe(ti,:);
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/norm(n);
        if normal*n' < 0
            ti = [ti(1),ti(3),ti(2)];
            n = -n;
        end
        
        switch ntri
            case 1
                color = 'b';
            case 2
                color = 'c';
            case 3
                color = 'r';
            case 4
                color = 'm';
        end
%         patch(Xe(ti,1),Xe(ti,2),Xe(ti,3),color,'EdgeColor','none')
%         quiver3(mean(Xe(ti,1)),mean(Xe(ti,2)),mean(Xe(ti,3)),n(1),n(2),n(3),0.1,'Color','b')

        % Add to surfh
        
        surfh(nTriEle+itri-1).iel  = iel;
        surfh(nTriEle+itri-1).Xe = Xe(ti,:);
        surfh(nTriEle+itri-1).xp = Xe(ti,1);
        surfh(nTriEle+itri-1).yp = Xe(ti,2);
        surfh(nTriEle+itri-1).zp = Xe(ti,3);
        surfh(nTriEle+itri-1).faceNormal = n;
        surfh(nTriEle+itri-1).ElementNormal = normal;

        % Assembly
        lo = (nTriEle+itri-1)*3-2;
        up = (nTriEle+itri-1)*3;
        surfX(lo:up,1:3) = Xe(ti,:);
        
    end

    nTriEle = nTriEle + ntri;

    end

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

        end

        % New order of polygon points
        k = convhull(Xe2D(:,1),Xe2D(:,2));
        k1 = k(1:end-1);
    end
