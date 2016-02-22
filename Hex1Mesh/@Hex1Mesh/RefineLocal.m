function [HangNodes, HangNodesM] = RefineLocal(T, ele)

if strcmpi(ele,'all')
    ele = 1:size(T.Connectivity,1);
end

%% Refine Hex Mesh Locally
% nodes = T.Connectivity;
% X = T.X;
xnod = T.XC;
ynod = T.YC;
znod = T.ZC;
% Number of nodes, keeps track of the latest node number
nmax = length(xnod);

nodes = T.Connectivity;
X = T.Points;

% Loop over all elements that are about to be refined
% For every element we will get 8 subelements
for iel = ele
    iv = T.Connectivity(iel,:); % local node numbers
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
    XN(XNeind,:) = (T.Points(edges(eind,1),:)+T.Points(edges(eind,2),:))/2;
    XN(XNMind,:) = [xm,ym,zm];
    XN(XNfind,:)=(T.Points(faces(faceind,1),:)+T.Points(faces(faceind,2),:)+T.Points(faces(faceind,3),:)+T.Points(faces(faceind,4),:))/4;
    
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
    INDS(J) = iv(I); % Existing Corner nodes
    
    
    NewNodes = (nmax+1:nmax+19)'; %Temporarly create
    
    %% Check if neighbors already have created same nodes
    % Temporarly add new coordinates to the list, new coordinates may
    % contain points that already exist
    X2 = [X;XN];
    % Those extra added points are found and stored as nodenumbers
    duplicateNodes = NewNodes( ismember(X2(NewNodes,:),X,'rows') );
    
    % The duplicates that we want to keep
    %     DupNodes = find(ismember(X,XN,'rows'));
    
    % First column stores indices to INDS second column stores node indices
    K = zeros(19,2);
    % We want to populate the noncorner nodes
    K(1:19,1) = [2,4,5,6,8,10,11,12,13,14,15,16,17,18,20,22,23,24,26]';
    % Loop over all 19 non corner nodes and look if the new coordinate
    % already exists in the domain. If it exists, add the index of the
    % existing node to the list, if it does not exist, increment the
    % nodenumber and add it to the list.
    for i = 1:19
        ind = find(ismember(X,XN(i,:),'rows'),1);
        if ~isempty(ind)
            K(i,2)=ind;
        else
            K(i,2)=nmax+1;
            nmax=nmax+1;
        end
    end
    
    % Fill the local INDS list
    INDS(K(:,1)) = K(:,2);
    
    % The new global nodes in a connectivity matrix, rady to be appended to
    % the rest of the Connectivity matrix
    C2 = INDS(locnodes);
    
    % Delete the duplicate points from the new Coordinate matrix.
    X2(duplicateNodes,:) = [];
    % Replace coord matrix
    X = X2;
    
    % Add new element matrix to the rest
    nodes = [nodes;C2];
    
    
    %% Viz mesh
    
    %     VizMeshTemp(T,nodes,X)
    
    %% Hanging nodes
    EdgeNodes = find(ismember(X,XN(XNeind,:),'rows'));
    FaceNodes = find(ismember(X,XN(XNfind,:),'rows'));
    MidNode = find(ismember(X,XN(XNMind,:),'rows'));
    
    T.Element(iel).HangNodes.EdgeNodes = EdgeNodes;
    T.Element(iel).HangNodes.FaceNodes = FaceNodes;
    T.Element(iel).HangNodes.MidNode = MidNode;
    T.Element(iel).HangNodes.ParentElement = iel;
    
    
    
end
% Remove the mother element, since it has been replaced by 8 smaller
% elements
nodes(ele,:) = [];

%% HangNodes handling
HangNodes = [T.Element.HangNodes];
ParentElement = [HangNodes.ParentElement];
EdgeNodes = [HangNodes.EdgeNodes];
FaceNodes = [HangNodes.FaceNodes];
MidNode = [HangNodes.MidNode];
HangNodesM = [ParentElement;EdgeNodes;FaceNodes;MidNode];
T.HangNodes = HangNodes;
% Hangnodes are ready to be returned
% Return HangNodes, HangNodesM

%% Update Face and edge matrices
Q = zeros(T.nele*6,4); %Stores faces
edges1 = zeros(T.nele*12,2); %stores edges

T.Element(T.nele).faces = []; %initialize T.Element.faces
T.Element(T.nele).edges = []; %initialize T.Element.edges

T.nele = size(nodes,1);
lof = 1; loe = 1;
for iel = 1:T.nele
    upf = lof+5;
    upe = loe+11;
    %Faces
    iface = [nodes(iel,[1,2,3,4]);...
        nodes(iel,[5,6,7,8]);...
        nodes(iel,[1,5,8,2]);...
        nodes(iel,[4,6,5,1]);...
        nodes(iel,[3,7,6,4]);...
        nodes(iel,[2,8,7,3])];
    Q(lof:upf,:) = iface;
    T.Element(iel).faces = iface;
    
    %edges
    n = nodes(iel,:);
    I = [1,2,4,3,5,8,6,7]; %Numbering order (indices)
    iedges = [n(I(1)),n(I(2));...
        n(I(2)),n(I(4));...
        n(I(4)),n(I(3));...
        n(I(3)),n(I(1));...
        n(I(1)),n(I(5));...
        n(I(3)),n(I(7));...
        n(I(4)),n(I(8));...
        n(I(2)),n(I(6));...
        n(I(6)),n(I(5));...
        n(I(5)),n(I(7));...
        n(I(7)),n(I(8));...
        n(I(8)),n(I(6))];
    edges1(loe:upe,:) = iedges;
    T.Element(iel).edges = iedges;
    %counter
    lof = upf +1;
    loe = upe +1;
end
T.Connectivity = nodes;
T.Points = X;
T.XC = X(:,1);
T.YC = X(:,2);
T.ZC = X(:,3);

T.Faces = Q;
T.nnod = length(T.XC);
T.edges = edges1;



end
