clear
close all
clc


x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
ne = 30;
nxe = ne;nye = ne;nze = ne;

disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

% hv = mesh1.vizMesh('ElementNumbers','NodeNumbers');

N1 = mesh1.Neighbors('Structured');


xnod = mesh1.XC;
ynod = mesh1.YC;
znod = mesh1.ZC;
%% Surface function
% Create the surface function
R = 0.89;
% surfaceFunction = @(x,y) ((x - x0).^2+(y - y0).^2).^.5-R;
% xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
xc = 0; yc = 0; zc = 0;
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;

phi = surfaceFunction(xnod, ynod, znod);


%% Surface
tic
[surfX, surfh, ExcessiveCurvatureInds] = mesh1.CutP1(phi);
CutP1Time = toc
h = mesh1.vizMesh([ExcessiveCurvatureInds.iel],78);
h.patch.EdgeColor = [0,0,1];
h.patch.FaceColor = [1,1,1];
% h.patch.FaceAlpha = 0.7;
mesh1.vizP1Surf()

%% Triangulation
XT = surfX;
isequalPoint = @(a,b) (sum(isequalAbs(a, b))==3);
%Use neighbor information to triangulate!
nTri = size(surfX,1)/3
tri = zeros(nTri,3);
% xfigure(1);
xfigure(1);
Q = zeros(nTri,1);
D(nTri).CommonNodes = 0;
AllHexEle = [surfh.iel]';
surfEle = mesh1.SurfaceInfo.CutElements;

tri = zeros(nTri, 3);
itri = 1;
xfigure; axis equal; hold on;
tic
tismember = [];
talternative = [];
for iP = 1:3:size(surfX,1)

    iel = surfh(itri).iel; % Mother element
    if itri == 1
         tri(1,1:3) = [1,2,3];
         xe = surfX(tri(itri,:),:);
%          patch(surfX(tri(itri,:),1),surfX(tri(itri,:),2),surfX(tri(itri,:),3),'m')
         
         %temp
         n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
         n = n/norm(n);
%          quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1)
         
         surfh(itri).triS = tri(itri,:);
         surfh(itri).faceNormal = n;
                    
         itri = itri +1;
    else
        
        %The overal strategy is to loop over each set of triangle points
        %and compare each point to the previous added list of points. If
        %the current point exists in the list, then it was previously added
        %and its position in the list is added to the list. Otherwise the
        %number of nodes is incremented and that number is added to the
        %triangulation list.
        % Example: 
        % First set is initiated to [1,2,3]
        % In the next iteration a point in the set is already in the list,
        % the point number is determined to be 2, so two new nodes are
        % added:
        % [1,2,3;
        %  4,5,3]
        % And so on...
        %
        %% Neighbors
        % Use neighborhood information to keep the search space small.
        
        %Get the current element neighbors and their neighbors to make sure
        %to capture all the neighbor nodes.
        neighHex = setdiff(unique([N1(N1(iel,N1(iel,:)>0),:);N1(iel,:)]),[0]);
%         t2=tic;
        neighHex = surfEle(ismember(surfEle,neighHex));    
%         tismember = [tismember;toc(t2)];
% %         bsxfun(@eq,surfEle,neighHex)
% %         length(ismember(surfEle,neighHex))
% %         length(sum(cell2mat(arrayfun(@(ind) surfEle == neighHex(ind), [1:length(neighHex)],'un',0)),2))
% %         isequal(ismember(surfEle,neighHex),sum(cell2mat(arrayfun(@(ind) surfEle == neighHex(ind), [1:length(neighHex)],'un',0)),2))
%         t1=tic;
%         neighHex= surfEle(logical(sum(cell2mat(arrayfun(@(ind) surfEle == neighHex(ind), [1:length(neighHex)],'un',0)),2)));
%         talternative = [talternative;toc(t1)];
        
        
        neighTri = find(ismember(AllHexEle,neighHex));
        neighPoints = sort([neighTri*3-2;neighTri*3-1;neighTri*3-0]);
        
%         mesh1.vizMesh(neighHex,3)
%         xx = reshape(surfX(neighPoints,1),3,[]);
%         yy = reshape(surfX(neighPoints,2),3,[]);
%         zz = reshape(surfX(neighPoints,3),3,[]);
%         patch(xx,yy,zz,'b')
        
%         xfigure(2)
        neighSpaceSize = length(neighPoints)-1;
        
        
        
        % iP is the first of the three points in the current set of
        % triangle points.
        % Loop over all three points in the set and compare each point to
        % all the previous points that were added. This set is named the
        % searchSpace.
        % It becomes obvious that as the loop progresses the searchSpace
        % increases and more comparisons are done, this makes this
        % algorithm scale bad. Or in other terms: the complexity is  ~O(n^1.5)
        % 
        % There is no need to compare with all the previous points. If the
        % neighborhood of the current set is known than the searchSpace is
        % kept pretty much constant. So we get to ~O(n) complexity
        
        % the local triangle node index is set to 1
        tind = 1;
        %% Loop over the three next points in the set of triangles
        
        for iPnt = iP:iP+2            
            % If the neighborhood point space is smaller than the space of
            % points we have found so far, then define the searchSpace as
            % the neighborhood points. Otherwise define the searchSpace as
            % the points we have found so far.
            % This is done because in the first few iterations the points
            % found so far are fewer than all the neighborhood points, but
            % later the neighborhgood is pretty mush constant where as the
            % points found so far would only increase.
            
%             if neighSpaceSize < iP+1
% %                 disp('neighborhood space')
% %                 neighPoints
%             else
% %                 disp('default space')
% %                 [1:iPnt-1]'
%             end
            PrevPoints = [1:iPnt-1]';
            SearchInds = PrevPoints(ismember(PrevPoints,neighPoints));
            
            searchSpace = surfX(SearchInds,:);
%             if neighSpaceSize < iP+1
%                 neighborHood = setdiff(neighPoints,iPnt)
%                 searchSpace = surfX(setdiff(neighPoints,iPnt),:);
%             else
% %                 [1:iPnt-1]'
%                 searchSpace = surfX(1:iPnt-1,:); 
%             end
            
            % Chose a point in the current set of trinagle points to measure the distance from
            Xp = surfX(iPnt,:);
            
            % Choose a set of points to measure the distance to. The set of
            % points is defined by the searchSpace.
            P = ones(size(searchSpace,1),1)*Xp;
            
            % The distance vectors between the chosen point and the searchSpace
            % points.
            D = searchSpace-P;
            
            % Eucledian distance
            NDist = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
            
            % Find the indices to the point that are identical (with room
            % for numerical error)
            indt = SearchInds(find(NDist<eps*100,1 ));
            
            % if duplicate points exist; set the local triangle index to
            % the index of the searchSpace and increment the triangle
            % index.
            % If no duplicate points are found; set the local triangle
            % index to the index of the current triangle set point
            if ~isempty(indt)
                tri(itri,tind) = indt;
                tind = tind +1;
            else
                tri(itri,tind) = iPnt;
                tind = tind +1;
            end
        end
        xe = surfX(tri(itri,:),:);
%         patch(surfX(tri(itri,:),1),surfX(tri(itri,:),2),surfX(tri(itri,:),3),'r')
        
        %temp
        n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
        n = n/norm(n);
%         quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1,'color','k')
        
        surfh(itri).triS = tri(itri,:);
        surfh(itri).faceNormal = n;
        
        
        itri = itri +1;
        
        
    end
    
end
tTriangulation = toc

%% Viz
xfigure(78); axis equal; hold on;
FV.Vertices = surfX;
FV.Faces = tri;
% shading interp
light
hp = patch(FV,'FaceColor','c','FaceLighting','gouraud');
view(136,16)
% hp.EdgeColor = 'none'











