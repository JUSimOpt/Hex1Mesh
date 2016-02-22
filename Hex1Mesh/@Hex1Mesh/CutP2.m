function [surfh,T] = CutP2(T, phi, level)

    nodes = T.Connectivity;
    xnod = T.XC;
    ynod = T.YC;
    znod = T.ZC;

    SurfEle = find(sum(phi(nodes)>0,2) > 0 & sum(phi(nodes)>0,2)< 8);
    nsurfEle = length(SurfEle);

    surfh(nsurfEle*4).iel = [];

    nTri = 1;
    % Loop over all Surface Elements
    for iel = SurfEle'
        %% Compute cut points
        iv = nodes(iel,:);

        %% Viz Faces
        %
        %     faces = T.Element(iel).faces;
        %     for j = 1:size(faces,1)
        %         ivv = faces(j,:);
        %         xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
        %         XC =[xc,yc,zc];
        %         XM = mean(XC,1);
        %         %         text(XM(1),XM(2),XM(3),num2str(j),'BackgroundColor','m')
        %     end

        edges = T.Element(iel).edges;
        %% Viz edges
        %     ed = unique(edges);
        %     for j = 1:8
        %         text(xnod(ed(j)),ynod(ed(j)),znod(ed(j)),num2str(ed(j)),'BackgroundColor','w')
        %     end

        %% Loop over all edges
        P = NaN(12,3);
        for j = 1:size(edges,1)
            ivv = edges(j,:);
            xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
            XC =[xc,yc,zc];
            q  = phi(edges(j,:));
            if sign(q(1)) ~= sign(q(2))
                [u,k] = sort(q);
                Xi = XC(k,:);
                P(j,:) = Xi(1,:)+(Xi(2,:)-Xi(1,:)).*((level-u(1))/(u(2)-u(1)));
            end
        end

        indP = find(~any(isnan(P),2));
        P = P(indP,:);
%         plot3(P(:,1),P(:,2),P(:,3),'ob')

        if size(P,1) == 3
            %% Tri
            np = 3;
            k1= RotatePointsToPlane(P);
            P = P(k1,:);
            E = edges(indP,:);
            E = E(k1,:);

            [Xe,tri,surfP2] = P2Points(T,iel,np,P,phi,iv,E,[xnod,ynod,znod],level);

            surfh(nTri).Xe = Xe;
            surfh(nTri).tri = surfP2.tri;
            surfh(nTri).iel = surfP2.iel;
            surfh(nTri).ElementNormal = surfP2.ElementNormal;
            nTri=nTri+1;
        elseif size(P,1) == 4
            %% Quad
            np = 4;
            k1= RotatePointsToPlane(P);
            P = P(k1,:);
            E = edges(indP,:);
            E = E(k1,:);   
            
            [Xe,tri,surfP2] = P2Points(T,iel,np,P,phi,iv,E,[xnod,ynod,znod],level);

            for i = 1:2
                surfh(nTri).Xe = Xe;
                surfh(nTri).tri = surfP2(i).tri;
                surfh(nTri).iel = surfP2(i).iel;
                surfh(nTri).ElementNormal = surfP2(i).ElementNormal;
                nTri=nTri+1;
            end

            

        elseif size(P,1) == 5
            np = 5;
            k1= RotatePointsToPlane(P);
            P = P(k1,:);
            E = edges(indP,:);
            E = E(k1,:);

            [Xe,tri,surfP2] = P2Points(T,iel,np,P,phi,iv,E,[xnod,ynod,znod],level);

            for i = 1:3
                surfh(nTri).Xe = Xe;
                surfh(nTri).tri = surfP2(i).tri;
                surfh(nTri).iel = surfP2(i).iel;
                surfh(nTri).ElementNormal = surfP2(i).ElementNormal;
                nTri=nTri+1;
            end

        elseif size(P,1) == 6
            np = 6;
            k1= RotatePointsToPlane(P);
            P = P(k1,:);
            E = edges(indP,:);
            E = E(k1,:);

            [Xe,tri,surfP2] = P2Points(T,iel,np,P,phi,iv,E,[xnod,ynod,znod],level);

            for i = 1:4
                surfh(nTri).Xe = Xe;
                surfh(nTri).tri = surfP2(i).tri;
                surfh(nTri).iel = surfP2(i).iel;
                surfh(nTri).ElementNormal = surfP2(i).ElementNormal;
                nTri=nTri+1;
            end
        else
            error('unknown element')
        end


    end
    
    T.SurfaceP2 = surfh;
    T.SurfaceP2Info.CutElements = SurfEle;
    T.SurfaceP2Info.NCutElements = length(SurfEle);

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

function [Xe, tri,surfP2] = P2Points(T,iel,np,P,phi,iv,E,X,level )
    %% Get mid points
    ed = [[1:np]',[2:np,1]'];
    Xe = zeros(np*2,3);
    xind = 1;
    L = zeros(np,1);
        
    for i=1:np
        xp = P(ed(i,:),:);
%         xfigure(45); axis equal; hold on;
%         plot3(xp(:,1),xp(:,2),xp(:,3),'-ob')
        Xe(xind,:)=xp(1,:);

        xind = xind+1;
        xm =  mean(xp,1);

        I = [1,2,4,3,5,8,6,7]; %Proper numbering
        u = phi(iv(I));
        
        %% In plane additional points
        
        %Face normal
        fn = cross(X(E(i,2),:)-X(E(i,1),:),xp(2,:)-xp(1,:));
        fn = fn/norm(fn);
        
        lambda = norm(xp(1,:)-xp(2,:))/2; %Length of the line

        xi = xm; %Start point for the search
        
        [fi, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
        PHI_xi = dot(u,fi);
%         text(xi(1),xi(2),xi(3),num2str(PHI_xi),'BackgroundColor','w')
        
        %search direction
        v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
        % In plane search direction
        pv = v - (dot(v,fn)/norm(fn)^2).*fn;
        
        npv = pv/norm(pv)*lambda;
%         quiver3(xi(1),xi(2),xi(3),npv(1),npv(2),npv(3),1,'Color','k')
%         quiver3(xi(1),xi(2),xi(3),-npv(1),-npv(2),-npv(3),1,'Color','k')

        xi1 = xi + npv;
        xi2 = xi - npv;
%         plot3(xi1(1),xi1(2),xi1(3),'xr')
%         plot3(xi2(1),xi2(2),xi2(3),'xr')
        
        [fi, ~, ~, ~, ~] = baseHexP1(T,iel,xi1);
        PHI_xi1 = dot(u,fi);
%         text(xi1(1),xi1(2),xi(3),num2str(PHI_xi1),'BackgroundColor','w')
        
        [fi, ~, ~, ~, ~] = baseHexP1(T,iel,xi2);
        PHI_xi2 = dot(u,fi);
%         text(xi2(1),xi2(2),xi2(3),num2str(PHI_xi2),'BackgroundColor','w')
        
        
        x1 = xi1;
        x2 = xi2;
        u1 = PHI_xi1;
        u2 = PHI_xi2;
        
        x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
%         [fi, fix, fiy, fiz, ~] = baseHexP1(T,iel,x3);
%         PHI_x3 = dot(u,fi);
        
        Xe(xind,:) = x3;
        xind = xind+1;
%         plot3(x3(1),x3(2),x3(3),'xr')
%         pause
        
    %     
    end
    
    
    
    xi = mean(Xe,1); %Start point for the search
    
    [fi, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
    v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
    nv = v/norm(v)*lambda;
    x1 = xi + nv;
    x2 = xi - nv;
    
    [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
    u1 = dot(u,fi);
    
    [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
    u2 = dot(u,fi);
    
    x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
    

    Xe(end+1,:) = x3;
%     plot3(x3(1),x3(2),x3(3),'xr')
%     pause


    E2 = E;
    [~,k] = sort(phi(E),2);
    for i=1:size(E,1)
        E2(i,:) = E(i,k(i,:));
    end
    en = mean(X(E2(:,2),:)-X(E2(:,1),:),1);
    en = en/norm(en);


    %% triangulation
    switch np
        case 3
            tri = [1,3,5,2,4,6];
            Xe(7,:)=[];
            surfP2.Xe = Xe;
            surfP2.tri = tri;
            surfP2.iel = iel;
            surfP2.ElementNormal = en;
        case 4
            tri = [1,3,5,2,4,9;1,5,7,9,6,8];
            
%             h = PlotP2(tri,Xe,20,'r',45);
%             for i = 1:size(Xe,1)
%                 text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
%             end
%             pause

            surfP2(1).tri = tri(1,:);
            surfP2(1).Xe = Xe(tri(1,:),:);
            surfP2(1).iel = iel;
            surfP2(1).ElementNormal = en;

            surfP2(2).tri = tri(2,:);
            surfP2(2).Xe = Xe(tri(2,:),:);
            surfP2(2).iel = iel;
            surfP2(2).ElementNormal = en;
        case 5
            %% Quad
            Xe(end,:)=[];


            xi = mean(Xe([1,5],:),1);
            lambda = norm(Xe(1,:)-Xe(5,:))/2;
            [~, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
            v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
            nv = v/norm(v)*lambda;
            x1 = xi + nv;
            x2 = xi - nv;
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
            u1 = dot(u,fi);
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
            u2 = dot(u,fi);
            x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
            Xe(end+1,:) = x3;

            xi = mean(Xe([1,7],:),1);
            lambda = norm(Xe(1,:)-Xe(7,:))/2;
            [~, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
            v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
            nv = v/norm(v)*lambda;
            x1 = xi + nv;
            x2 = xi - nv;
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
            u1 = dot(u,fi);
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
            u2 = dot(u,fi);
            x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
            Xe(end+1,:) = x3;

            tri = [1,3,5,2,4,11;1,5,7,11,6,12;1,7,9,12,8,10];

    %         h = TriP2(tri,Xe,20,'r',45);
    %         for i = 1:size(Xe,1)
    %             text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
    %         end
    %         return

            for i=1:3
                surfP2(i).tri = tri(i,:);
                surfP2(i).Xe = Xe(tri(i,:),:);
                surfP2(i).iel = iel;
                surfP2(i).ElementNormal = en;
            end

        case 6
            Xe(end,:)=[];

            xi = mean(Xe([3,7],:),1);
            lambda = norm(Xe(3,:)-Xe(7,:))/2;
            [~, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
            v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
            nv = v/norm(v)*lambda;
            x1 = xi + nv;
            x2 = xi - nv;
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
            u1 = dot(u,fi);
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
            u2 = dot(u,fi);
            x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
            Xe(end+1,:) = x3;

            xi = mean(Xe([1,9],:),1);
            lambda = norm(Xe(1,:)-Xe(9,:))/2;
            [~, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
            v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
            nv = v/norm(v)*lambda;
            x1 = xi + nv;
            x2 = xi - nv;
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
            u1 = dot(u,fi);
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
            u2 = dot(u,fi);
            x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
            Xe(end+1,:) = x3;

            xi = mean(Xe([1,3,7,9],:),1);
            [~, fix, fiy, fiz, ~] = baseHexP1(T,iel,xi);
            v = [dot(u,fix), dot(u,fiy), dot(u,fiz)];
            nv = v/norm(v)*lambda;
            x1 = xi + nv;
            x2 = xi - nv;
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x1);
            u1 = dot(u,fi);
            [fi, ~, ~, ~, ~] = baseHexP1(T,iel,x2);
            u2 = dot(u,fi);
            x3 = x1+(x2-x1).*((level-u1)/(u2-u1));
            Xe(end+1,:) = x3;

    %         for i = 1:size(Xe,1)
    %             text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
    %         end

            tri = [3,5,7,4,6,13; 3,7,9,13,8,15; 3,9,1,15,14,2; 1,9,11,14,10,12];


            for i=1:4
                surfP2(i).tri = tri(i,:);
                surfP2(i).Xe = Xe(tri(i,:),:);
                surfP2(i).iel = iel;
                surfP2(i).ElementNormal = en;
            end


            return
        otherwise
            error('unknown element!')
    end

end

