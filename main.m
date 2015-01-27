clear
close all
clc


x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;
nxe = 1;
nye = 1;
nze = 1;

disp('creating mesh...')
tic
mesh1 = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
toc

hv = mesh1.vizMesh('ElementNumbers','NodeNumbers');


locnodes = mesh1.Connectivity;
xnod = mesh1.xnod; ynod = mesh1.ynod; znod = mesh1.znod;
X= [xnod,ynod,znod];
% mh = mesh1.mesh


N1 = mesh1.Neighbors('Structured');


% N2 = mesh1.Neighbors('Naive');


X = [mesh1.xnod,mesh1.ynod,mesh1.znod];


%%



%% Refine Hex Mesh Locally
T = mesh1;
% set(hv.EleText,'Visible','off')
nele = [1,3]

% Number of nodes, keeps track of the latest node number
nmax = length(xnod);
% Initialize Hanging nodes matrix
HangNodes = [];
nodes = T.Connectivity;
X = T.X;
xnod = T.xnod;
ynod = T.ynod;
znod = T.znod;

% Loop over all elements that are about to be refined
% For every element we will get 8 subelements
for iel = nele
   
   iv = nodes(iel,:); % local node numbers
   xc = xnod(iv); % local coordinates
   yc = ynod(iv);
   zc = znod(iv);
   
   xm = mean(xc);ym = mean(yc);zm = mean(zc); %mid points of element
   
   edges = T.Element(iel).edges; % Edge list of current element.
   faces = T.Element(iel).faces; % Face list of current element. 
   %% Create New Coords
   % We start by creating a list of new points. The list includes points on
   % the mid points of edges, faces and element. 19 points total
   
   % Edge indices, the order is chosen to create the new points in the same
   % order as the mother nodes are numbered (along positive y, pos x and
   % then pos z. See http://www.mirzacenanovic.com/wp-content/uploads/2015/01/2015-01-09-14_50_26-Figure-1_-Hex1Mesh.png
   eind = [1,4,2,3,5,8,6,7,9,10,12,11];
   faceind = [1,3,4,6,5,2];
   %The order in which the subelement points are created is the same order
   %as for the mother element
   XNeind = [1,2,4,5,6,8,12,14,15,16,18,19]; %edge inds
   XNMind = 10; %mid
   XNfind = [3,7,9,11,13,17]; %face
   
   % Notice that we're creating 19 points regardless if these points exists
   % in neighboring elements. We'll deal with that later.
   XN = zeros(19,3);   % 19 New points
   XN(XNeind,:) = (X(edges(eind,1),:)+X(edges(eind,2),:))/2;
   XN(XNMind,:) = [xm,ym,zm];
   XN(XNfind,:)=(X(faces(faceind,1),:)+X(faces(faceind,2),:)+X(faces(faceind,3),:)+X(faces(faceind,4),:))/4;
    
   %Local node numbering same as for mother elements
   locnodes = [1  2  5  4  10 13 14 11;...
               2  3  6  5  11 14 15 12;...
               4  5  8  7  13 16 17 14;...
               5  6  9  8  14 17 18 15;...
               10 11 14 13 19 22 23 20;...
               11 12 15 14 20 23 24 21;...
               13 14 17 16 22 25 26 23;...
               14 15 18 17 23 26 27 24];
   %% Newinds  
   I = [1,2,4,3,5,8,6,7];
   J = [1,3,7,9,19,21,25,27];
   INDS = zeros(27,1); %Indices
   INDS(J) = iv(I); % Exiosting Corner nodes
   NewNodes = (nmax+1:nmax+19)'; %Temporarly create
   
   %% Check if neighbors already have created same nodes
   % Temporarly add new coordinates to the list, new coordinates may
   % contain points that already exist
   X2 = [X;XN]; 
   % Those extra added points are found and stored as nodenumbers
   duplicateNodes = NewNodes( ismember(X2(NewNodes,:),X,'rows') );
   % The duplicates that we want to keep
   DupNodes = find(ismember(X,XN,'rows')); 
   % These are all the nodes we want to keep, including the mother element
   % nodes
   PrescNodes = [iv(I)'; DupNodes];
   % We find out the new number of nodes to add
   nNewNodes = 27-length([iv(I)'; DupNodes]);
   % List of new nodes
   NewUniqueNodes = (nmax+1:nmax+nNewNodes)';
   % Update the max node number
   nmax = nmax+nNewNodes;
   % Update list of locally refined nodes with the new nodes, including the
   % duplicate ones
   INDS(setdiff(1:27,J)) = NewNodes;
   % The above step is done so that we can find out the indices of the
   % duplicate nodes
   DupNodeInds = find(ismember(INDS,duplicateNodes));
   % The prescribed node inds, includes the mother element nodes
   presc = [J';DupNodeInds];
   % Now we reset INDS and updated it with the New Unique Node numbers
   INDS = zeros(27,1); %Indices
   % We keep the corner nodes and the neighbor nodes
   INDS(presc) = [iv(I)'; DupNodes];
   % And add the new Unique Nodes in the rest of the positions
   free = setdiff(1:27,presc)';
   INDS(free) = NewUniqueNodes;
   
   % The new global nodes in a connectivity matrix, rady to be appended to
   % the rest of the Connectivity matrix
   C2 = INDS(locnodes);
   % Delete the duplicate points from the new Coordinate matrix.
   X2(duplicateNodes,:) = [];
   % Replace coord matrix
   X = X2;
   % Remove the mother element, since it has been replaced by 8 smaller
   % elements
   nodes(iel,:) = [];
   % Add new element matrix to the rest
   nodes = [nodes;C2];
   
   
   
        
   
   %% Hanging nodes
   
%    DupNodes = find(ismember(X,XN,'rows')) %Coordinates to Keep
   
   
   
%    duplicateNodes = NewNodes(find(ismember(X2(NewNodes,:),X,'rows'))) %Coordinates to replace with DupNodes
%    INDS(ismember(INDS,duplicateNodes)) = DupNodes;
   
   
   
   
   
%    unique(C2);
%    
%    numberOfUniqueNewNodes = sum(~ismember(X2(unique(C2),:), X,'rows'))
%    
%    size(NewNodes)
   
%    
   
   
   
   
   
%    DupNodes = setdiff(HangNodes,(1:length(X))') 
%    X2(DupNodes,:) = [];
%    HangNodes = setdiff(NewNodes,HangNodes); 
   
   
   
   
   

   
   
%    for j=1:size(X,1)
%        text(X(j,1),X(j,2),X(j,3),num2str(j),'BackgroundColor','w')
%    end
%    
   %% Viz mesh

   
   
end


% nodes
% X(47,:)=[];
% size(X)
% for i = 47:max(nodes(:))
%     nodes(nodes==i)=i-1;
% end
nodes
size(X)
max(nodes(:))

hv2 = mesh1.vizMesh();
for i = 1:size(nodes,1)
    iv = nodes(i,:);
    xc = X(nodes(i,:),1);
    yc = X(nodes(i,:),2);
    zc = X(nodes(i,:),3);
    
    for j=1:8
        text(xc(j),yc(j),zc(j),num2str(iv(j)),'BackgroundColor','w')
    end
end

