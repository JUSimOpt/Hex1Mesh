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
        XC
        YC
        ZC
        nnod
        nele
        edges
        Neighs
        Element
        HangNodes
        SurfaceP1
        SurfaceP1Points
        SurfaceP1Triangulation
        SurfaceP1Info
        
        SurfaceP2
        SurfaceP2Points
        SurfaceP2Triangulation
        SurfaceP2Info
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
            T.XC = T.X(:,1); T.YC = T.X(:,2); T.ZC = T.X(:,3);
            
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
            
            T.nnod = length(T.XC);
            T.edges = edges1;
            
            T.HangNodes = [];
            T.Element(1).HangNodes = [];
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
end

