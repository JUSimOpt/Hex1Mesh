clear
close all
clc


%% Mesh
x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 4;
nye = 4;
nze = 4;


disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc
hv = mesh1.vizMesh();

neighs = mesh1.Neighbors('Structured');

xnod = mesh1.xnod;
ynod = mesh1.ynod;
znod = mesh1.znod;
nodes = mesh1.Connectivity;
%     8-----7
%    /|    /|
%   6-----5 |
%   | 4...|.3
%   |/    |/
%   2-----1

%% Surface function
% Create the surface function
R = 0.89;
% surfaceFunction = @(x,y) ((x - x0).^2+(y - y0).^2).^.5-R;
% xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
xc = 0; yc = 0; zc = 0;
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;

%% Create surfaceClass
phi = surfaceFunction(xnod, ynod, znod);
SurfEle = find(sum(phi(nodes)>0,2) > 0 & sum(phi(nodes)>0,2)< 8);

%% Get cut elements
tic
% Loop over all elements
% Loop over all faces
%Loop over all edges
%Identify which face is cut
% what case

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
xfigure(2); hold on
light
for iel = SurfEle'
    iv = nodes(iel,:);
    
    
    faces = mesh1.mesh(iel).faces;
    
    Pe = NaN(2*6,3);
    ed = 1;
    normal = zeros(12,3);
    for j = 1:size(faces,1)
        ivv = faces(j,:);
        xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
        XC =[xc,yc,zc];
        q  = phi(faces(j,:));
        patch(xc,yc,zc,'w','FaceColor','none');  hold on;view(23,37)
%         for in = 1:4
%             text(xc(in),yc(in),zc(in), [num2str(ivv(in)),', ',num2str(q(in))] )
%         end
        axis equal tight;
        
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
    normal = sum(normal,1);
    normal = normal/norm(normal);
    
    

    
    Xe = unique(Pe,'rows','stable');
    Xe = Xe(~any(isnan(Xe),2),:);
    if size(Xe,1) == 4
        %       disp('quad')
        %       plot3(Xe(:,1),Xe(:,2),Xe(:,3),'b*')
        xp1 = Pe(1,:);
        xp2 = Pe(2,:);
        for i = 3:12
            xpi = Pe(i,:);
            if isequalPoint(xpi,xp2)
                if isodd(i)
                    xp3 = Pe(i+1,:);
                else
                    xp3 = Pe(i-1,:);
                end
            elseif isequalPoint(xpi,xp1)
                if isodd(i)
                    xp4 = Pe(i+1,:);
                else
                    xp4 = Pe(i-1,:);
                end
            end
        end
        
        Xe = [xp1;xp2;xp3;xp4];
        for i = 1:4
            text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i))
        end
        
%         quiver3(mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3)),normal(1),normal(2),normal(3),0.1)
        
        % Quad to triangles
        Xe1 = Xe(1:3,:);
        Xe2 = [Xe(3,:);Xe(4,:);Xe(1,:)];
        
        patch(Xe1(:,1),Xe1(:,2),Xe1(:,3),'w')
        patch(Xe2(:,1),Xe2(:,2),Xe2(:,3),'w')
        
        
        e1 = Xe(2,:)-Xe(1,:);
        e2 = Xe(3,:)-Xe(1,:);
        e3 = Xe(4,:)-Xe(1,:);
        fN1 = cross(e2,e1);
        fN1 = fN1/norm(fN1); 
        fN2 = cross(e3,e1);
        fN2 = fN2/norm(fN2); 
        
        quiver3(mean(Xe1(:,1)),mean(Xe1(:,2)),mean(Xe1(:,3)),fN1(1),fN1(2),fN1(3),0.1)
        quiver3(mean(Xe2(:,1)),mean(Xe2(:,2)),mean(Xe2(:,3)),fN2(1),fN2(2),fN2(3),0.1)
        
        
        %         Add to surfh
        surfh(nTriEle).iel  = iel;
        surfh(nTriEle).Xe    = Xe1;
        surfh(nTriEle).xp = Xe1(:,1);
        surfh(nTriEle).yp = Xe1(:,2);
        surfh(nTriEle).zp = Xe1(:,3);
        surfh(nTriEle).nei = neighs(iel,:);
        surfh(nTriEle).faceNormal = fN1;
        
        surfh(nTriEle + 1).iel  = iel;
        surfh(nTriEle + 1).Xe    = Xe2;
        surfh(nTriEle + 1).xp = Xe2(:,1);
        surfh(nTriEle + 1).yp = Xe2(:,2);
        surfh(nTriEle + 1).zp = Xe2(:,3);
        surfh(nTriEle + 1).nei = neighs(iel,:);
        surfh(nTriEle + 1).faceNormal = fN2;
        
        %Assembly
        lo = nTriEle*3-2;
        up = nTriEle*3;
        surfX(lo:up,1:3) = Xe1;
        
        lo = (nTriEle+1)*3-2;
        up = (nTriEle+1)*3;
        surfX(lo:up,1:3) = Xe2;
        
        nTriEle = nTriEle +2;
        

    elseif size(Xe,1) == 3
        %        disp('tri')
        %        plot3(Xe(:,1),Xe(:,2),Xe(:,3),'b*')
    end
    
    pause
    
    
    
    
end
axis equal tight
toc