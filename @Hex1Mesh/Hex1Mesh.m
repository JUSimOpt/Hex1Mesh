classdef Hex1Mesh < matlab.mixin.Copyable
    %Hex1Mesh Linear structured Hexahedral mesh
    % 
    % Hexahedral is numbered bottom face starting in lower right corner going
    % counter-clockwise viewing from the bottom. Bottom to top. 
    %
    %   8-----7
    %  /|    /|
    % 6-----5 |
    % | 4...|.3
    % |/    |/ 
    % 2-----1
    %
    % Faces are numbered in the following way
    % face 1 = 1,2,4,3
    % face 2 = 5,7,8,6
    % face 3 = 1,5,6,2
    % face 4 = 3,7,5,1
    % face 5 = 4,8,7,3
    % face 6 = 2,6,8,3

    properties
        Connectivity
        Faces
        xnod
        ynod
        znod
        nnod
        nele
        edges
        Neighs
        Element
        HangNodes
    end
    
    properties (Hidden)
        X
        nx
        ny
        nz
    end
    
    properties (Access = private)
        x0
        x1
        y0
        y1
        z0
        z1
        nxe
        nye
        nze
    end
    
    methods
        function T = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
            
            nele = nxe*nye*nze;
            if nele > 100000
                btn = questdlg(['Warning, ',num2str(nele),' elements are about to be created. Do you wish to continue?'], ...
                    '¡Ay! ¡Mucho elementos!!','Si','No, no, nooo', 'No, no, nooo');
                putvar(btn)
                
                switch btn
                    case 'Si'
                        disp('Creating huge amounts of elements...')
                    case 'No, no, nooo'
                        disp('No mesh is created!')
                        return
                end
            end
            
            x = linspace(x0,x1,nxe+1);
            y = linspace(y0,y1,nye+1);
            z = linspace(z0,z1,nze+1);
            [MX, MY, MZ] = meshgrid(x,y,z);
            
            T.X = [MX(:),MY(:),MZ(:)];
            
            nx = nxe+1;
            ny = nye+1;
            nz = nze+1;
            
            T.nxe = nxe;
            T.nye = nye;
            T.nze = nze;
            
            T.nele = nxe*nye*nze;
            
            nodes = zeros(T.nele, 8);
            c = 1;
            iel = 1;
            for k = 1:nz
                for j = 1:nx
                    for i = 1:ny
                        
                        if i ~= ny && j ~=nx && k ~= nz
                            nodes(iel,:) = [c,c+1,c+ny+1,c+ny,c+nx*ny,c+nx*ny+ny,c+nx*ny+ny+1,c+nx*ny+1];
                            iel = iel + 1;
                        end
                        c = c+1;
                    end
                end
            end
            
            T.nx = nx;
            T.ny = ny;
            T.nz = nz;
            
            
            T.Connectivity = nodes;
            T.xnod = T.X(:,1); T.ynod = T.X(:,2); T.znod = T.X(:,3);
            
            Q = zeros(T.nele*6,4);
            edges1 = zeros(T.nele*12,2);
            
            T.Element(T.nele).faces = [];
            T.Element(T.nele).edges = [];
            
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
                I = [1     2     4     3     5     8     6     7];
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
            
            T.Faces = Q;
            
            T.nnod = length(T.xnod);
            T.edges = edges1;
            
            T.HangNodes = [];
            T.Element(1).HangNodes = [];
        end
        
        function [neighs,T] = Neighbors(T,varargin)
            
            if nargin > 1
                searchmethod = varargin{1};
            else
                searchmethod = 'Structured'
            end
            
            switch searchmethod
                case 'Structured'
                    %% Structured
                    neighs = zeros(T.nele,6);
                    iel = 1;
                    nx = T.nxe;
                    ny = T.nye;
                    nz = T.nze;
                    for k = 1:T.nze
                        for j = 1:nx
                            for i = 1:T.nye
                                
                                %y-direction
                                if i == ny && i == 1
                                    n1 = 0;
                                    n2 = 0;
                                elseif i == ny
                                    n1 = 0;
                                    n2 = iel-1;
                                elseif i == 1
                                    n2 = 0;
                                    n1 = iel+1;
                                else
                                    n1 = iel+1;
                                    n2 = iel-1;
                                end
                                
                                %x-direction
                                if j==nx && j==1
                                    n3 = 0;
                                    n4 = 0;
                                elseif j==nx
                                    n3 = 0;
                                    n4 = iel-ny;
                                elseif j == 1
                                    n3 = iel+ny;
                                    n4 = 0;
                                else
                                    n3 = iel+ny;
                                    n4 = iel-ny;
                                end
                                
                                %z-direction
                                if k==nz && k==1
                                    n5 = 0;
                                    n6 = 0;
                                elseif k == nz
                                    n5 = 0;
                                    n6 = iel-ny*nx;
                                elseif k == 1
                                    n5 = iel+ny*nx;
                                    n6 = 0;
                                else
                                    n5 = iel+ny*nx;
                                    n6 = iel-ny*nx;
                                end
                                
                                neighs(iel,:) = [n1,n2,n3,n4,n5,n6];
                                T.Element(iel).neighs = [n1,n2,n3,n4,n5,n6];
                                
                                iel = iel +1;
                                
                            end
                        end
                    end
                    
                case 'Naive'
                    %
                    %   8-----7
                    %  /|    /|
                    % 6-----5 |
                    % | 4...|.3
                    % |/    |/
                    % 2-----1
                    nodes = T.Connectivity;
                    X = T.X;
                    
                    
                    nele = size(nodes,1);
                    xnod = X(:,1); ynod = X(:,2); znod = X(:,3);
                    
                    nnod = size(xnod,1);
                    
                    %% Element size
                    EleVolume = zeros(nele,1);
                    EleSize = zeros(nele,3);
                    for iel = 1:nele
                        iv = nodes(iel,:);
                        xc = X(iv,1); yc = X(iv,2); zc = X(iv,3);
                        
                        %     EleVolume(iel) =  HexVolume(iv,X);
                        
                        xm = mean(X(iv,1)); ym = mean(X(iv,2)); zm = mean(X(iv,3));
                        
                        dxm = max(abs(xm-xc));
                        dym = max(abs(ym-yc));
                        dzm = max(abs(zm-zc));
                        EleSize(iel,:) = [dxm,dym,dzm];
                    end
                    ms = 2*norm([max(EleSize(:,1)),max(EleSize(:,2)),max(EleSize(:,3))]);
                    
                    %% Find neighbors
                    neighs = zeros(nele,6);
                    ElementSpace = 1:nele;
                    for iel = 1:nele
                        iv = nodes(iel,:);
                        xc = X(iv,1); yc = X(iv,2); zc = X(iv,3);
                        
                        %% Hex faces
                        f1 = nodes(iel,[1,2,3,4]);
                        f2 = nodes(iel,[5,6,7,8]);
                        f3 = nodes(iel,[1,5,8,2]);
                        f4 = nodes(iel,[4,6,5,1]);
                        f5 = nodes(iel,[3,7,6,4]);
                        f6 = nodes(iel,[2,8,7,3]);
                        
                        %% Mid point of element iel
                        
                        xm = mean(X(iv,1));
                        ym = mean(X(iv,2));
                        zm = mean(X(iv,3));
                        
                        %% Loop over all neighbor elements
                        n1 = 0;
                        n2 = n1; n3=n1;n4=n1;n5=n1;n6=n1;
                        n = [0,0,0,0,0,0];
                        for j = ElementSpace
                            if j == iel
                                continue
                            end
                            ivj = nodes(j,:);
                            xmj = mean(X(ivj,1));
                            ymj = mean(X(ivj,2));
                            zmj = mean(X(ivj,3));
                            dijx = abs(xmj-xm);
                            dijy = abs(ymj-ym);
                            dijz = abs(zmj-zm);
                            
                            dij = sqrt((xmj-xm)^2+(ymj-ym)^2+(zmj-zm)^2);
                            
                            if dij > ms
                                continue
                            end
                            %        j
                            %% Get faces of current neighbordhood element
                            fj = [nodes(j,[1,2,3,4]);...
                                nodes(j,[5,6,7,8]);...
                                nodes(j,[1,5,8,2]);...
                                nodes(j,[4,6,5,1]);...
                                nodes(j,[3,7,6,4]);...
                                nodes(j,[2,8,7,3])];
                            
                            %        [f1;f2;f3;f4;f5;f6]
                            %        fj
                            ind = find(sum(ismember([f1;f2;f3;f4;f5;f6], fj),2) == 4);
                            if ~isempty(ind)
                                n(ind) = j;
                            end
                            
                            
                        end
                        %     iel
                        %     n
                        %     pause
                        neighs(iel,:) = n;
                        
                    end
            end
            
            T.Neighs = neighs;
        end
        
        function h = vizMesh(T,varargin)
            % h = vizMesh()
            % h = vizMesh(ele,properties)
            % ele is a list of elements
            % properties:
            % 'NodeNumbers'
            %
            ele = 1:size(T.Connectivity,1);
            if nargin > 1
                if isa(varargin{1},'double')
                    ele = varargin{1};
                end
            end
            
            h.fig = xfigure;
            
            ele = ele(:);
            fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
            
            h.patch = patch(T.xnod(T.Faces(fele(:),:)'),T.ynod(T.Faces(fele(:),:)'),T.znod(T.Faces(fele(:),:)'),'w','FaceColor','none');
            xlabel('X'); ylabel('Y'); zlabel('Z')
            axis equal tight
            set(h.fig,'name','Hex1Mesh')
            view(-55,45)
            
            if isenabled('NodeNumbers',varargin)
                if T.nele > 100
                    warning('Cannot draw NodeNumbers, too many elements')
                    return
                end
                h.NodeText = [];     
%                 unique(T.Connectivity(:),'stable')'
%                 1:T.nnod
                for i = 1:T.nnod
                    h.NodeText = [h.NodeText; text(T.xnod(i),T.ynod(i),T.znod(i),num2str(i),'BackgroundColor','w') ];
                end
            end
            
            if isenabled('ElementNumbers',varargin)
                if T.nele > 100
                    warning('Cannot draw ElementNumbers, too many elements')
                    return
                end
                
                h.EleText = [];
                for i = sort(ele)'
                    xm = mean(T.xnod(T.Connectivity(i,:)));
                    ym = mean(T.ynod(T.Connectivity(i,:)));
                    zm = mean(T.znod(T.Connectivity(i,:)));
                    h.EleText = [h.EleText; text(xm,ym,zm,num2str(i),'BackgroundColor','y')];
                end
                
            end
            
            
            
        end
        
        function V = ElementVolume(T,iel)
            % Efficient Computation of Volume of Hexhedral Cells by J.Grandy
            node = T.Connectivity(iel,:);
            X = T.X;
            v11 = X(node(6),:)-X(node(2),:);
            v12 = X(node(1),:)-X(node(2),:);
            v13 = X(node(4),:)-X(node(5),:);
            
            v21 = v11;
            v22 = X(node(8),:)-X(node(2),:);
            v23 = X(node(5),:)-X(node(7),:);
            
            v31 = v11;
            v32 = X(node(3),:)-X(node(2),:);
            v33 = X(node(7),:)-X(node(4),:);
            
            V = 1/6*(det([v11',v12',v13'])+det([v21',v22',v23'])+det([v31',v32',v33']));
        end
        
        
        
    end
    
    methods (Access = private)
        
        function rv = isenabled(mode, varargin)
            if nargin < 1
                error('No arguments')
            end
            varargin = varargin{:};
            
            ind = find(strcmpi(varargin,mode), 1);
            if ~isempty(ind)
                rv = 1;
            else
                rv = 0;
            end
        end
        
    end
    
    %Hide some of the inherited methods from handle
    methods(Hidden)
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function delete(varargin)
            delete@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
    end
    
end

