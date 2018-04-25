function [tri,X,T] = TriangulateP1(T)
%TRIANGULATEP1 Creates a triangulation from a set of P1 triangles
%   Triangulates the P1 Surface triangle set surfX
%   
%   [tri,X] = TriangulateP1()
%   tri is the triangulation matrix m-by-3
%   X is the set of triangle coordinates same as (surfX)
%   Usage:
%   T.TriangulateP1( ... )
%   or
%   TriangulateP1(T, ... )
%   T is the Hex1Mesh class

surfX = T.SurfaceP1Points;
surfh = T.SurfaceP1;

AllHexEle = [surfh.iel]';
surfEle = T.SurfaceP1Info.CutElements;
if isempty(T.Neighs)
    N1 = T.Neighbors('Structured');
end
N1 = T.Neighs;

nTri = size(surfX,1)/3;
tri = zeros(nTri, 3);
itri = 1;
for iP = 1:3:size(surfX,1)

    iel = surfh(itri).iel; % Mother element
    if itri == 1
         tri(1,1:3) = [1,2,3];
         xe = surfX(tri(itri,:),:);
         
%          xfigure;
%          patch(xe(:,1),xe(:,2),xe(:,3),'m'); hold on; axis equal; light
         
         %temp
         n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
         n = n/norm(n);
%          quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1)
         
         surfh(itri).triS = tri(itri,:);
         surfh(itri).faceNormal = n;
                    
         itri = itri +1;
%          pause
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

        neighHex = surfEle(ismember(surfEle,neighHex));    

        
        
        neighTri = find(ismember(AllHexEle,neighHex));
        neighPoints = sort([neighTri*3-2;neighTri*3-1;neighTri*3-0]);
        
%         neighSpaceSize = length(neighPoints)-1;
        
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
            
            PrevPoints = [1:iPnt-1]';
            SearchInds = PrevPoints(ismember(PrevPoints,neighPoints));
            
            searchSpace = surfX(SearchInds,:);
            
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
        
%         patch(xe(:,1),xe(:,2),xe(:,3),'m'); axis equal;
        

        n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
        n = n/norm(n);
        
%         quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1)
        
        surfh(itri).triS = tri(itri,:);
        surfh(itri).faceNormal = n;
        
        
        itri = itri +1;
%         pause
        
    end
    
end

%% Weld
X = surfX;

mt = 0;
map = zeros(size(X,1),2);
c = 0;
for iel = 1:size(tri)
    itri = tri(iel,:);
    
    if iel == 1
        for it = itri
            mt=mt+1;
        end
        map(1:3,1)=[1:3]';
        map(1:3,2)=[1:3]';
        c = 3;
        mt = 4;
    else
        for indt=1:3
            it = itri(indt);
            if it >= mt && ~any(it==map(:,1))
                map(c+1,:) = [it,mt];
                mt=mt+1;
                c=c+1;
            else
            end
            
        end
        
    end
end
map = map(all(map~=0,2),:);

for im = 1:size(map,1)
    tri(tri==map(im,1)) = map(im,2);
end
X =  X(map(:,1),:);

%% Return
T.SurfaceP1Triangulation = tri;
T.SurfaceP1 = surfh;
end
