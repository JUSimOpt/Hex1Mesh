function  VizMeshTemp(T,nodes,X)

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
T.X = X;
T.xnod = X(:,1);
T.ynod = X(:,2);
T.znod = X(:,3);

T.Faces = Q;
T.nnod = length(T.xnod);
T.edges = edges1;


%% Vix mesh
% nodes
% size(X)
% max(nodes(:))
hv2 = T.vizMesh();
for i = 1:size(nodes,1)
    iv = nodes(i,:);
    xc = X(nodes(i,:),1);
    yc = X(nodes(i,:),2);
    zc = X(nodes(i,:),3);
    
    for j=1:8
        text(xc(j),yc(j),zc(j),num2str(iv(j)),'BackgroundColor','w')
    end
end

end

