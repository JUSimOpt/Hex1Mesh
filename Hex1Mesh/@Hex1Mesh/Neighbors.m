function [neighs,T] = Neighbors(T,varargin)
%Neighbors Computes neighbors for the Hex1Mesh
%   [neighs,T] = Neighbors(T)
%   [neighs,T] = Neighbors(T,searchmethod)
%   Generates an m-by-6 matrix containing indices to neighbor element numbers for
%   each side of the element - stored in neighs.
%   T is the Hex1Mesh class
%   neighs is the returened neighbor matrix
%   searchmethod is a string either 'Structured' or 'Naive'. Default is
%   Structured and there is no reason to use Naive
%
%   Usage:
%   T.Neighbors( ... )
%   or
%   Neighbors(T, ... )
%   T is the Hex1Mesh class


if nargin > 1
    searchmethod = varargin{1};
else
    searchmethod = 'Structured';
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