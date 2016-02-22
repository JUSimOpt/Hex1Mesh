function [tri,surfX,T] = TriangulateP2(T)
%TRIANGULATEP2 Creates a triangulation from a set of P2 triangles
%   Triangulates the P2 Surface triangle set surfX
%   
%   [tri,X] = TriangulateP2()
%   tri is the triangulation matrix m-by-6
%   X is the set of triangle coordinates same as (surfX)
%   Usage:
%   T.TriangulateP2( ... )
%   or
%   TriangulateP2(T, ... )
%   T is the Hex1Mesh class

    if isempty(T.SurfaceP2)
        error('Surface must exist!')
    end
    
    if isempty(T.SurfaceP2Triangulation)
        surfh = T.SurfaceP2;
        X = cell2mat({surfh.Xe}');
        tri = cell2mat({surfh.tri}');
        tri = tri*0;
        X = X*0;

        % tri=[];
        % X=[];
        % iel0= 1;
        itri = 1;
        xlo = 1;
        % maxt = 0;
        for i=1:length(surfh)
            %     i
            Xe = surfh(i).Xe;
            tt = surfh(i).tri;
            iel1 = surfh(i).iel;

            if isempty(tt)
                break
            end

            %     ntri = size(tt,1);

            if i==1
                iel0 = iel1;
                %         X = [X;Xe];
                xup = xlo+size(Xe,1);
                X(xlo:xup-1,:)=Xe;
                xlo = xup;
                maxt =0;
            end

            if iel1 ~= iel0
                %Unique element
                %         X = [X;Xe];
                xup = xlo+size(Xe,1);
                X(xlo:xup-1,:)=Xe;
                xlo = xup;

                iel0 = iel1;

                maxt = max(tri(:));
                %         tri = [tri;tt+maxt];
                tri(itri,:) = tt+maxt;
            else
                %         tri = [tri;tt+maxt];
                tri(itri,:) = tt+maxt;
            end
            itri = itri+1;

            % pause


        end
        % max(tri(:))
        
        X = X(1:max(tri(:)),:);
        
        T.SurfaceP2Triangulation = tri;
        T.SurfaceP2Points = X;
    end

    tri = T.SurfaceP2Triangulation ;
    surfX = T.SurfaceP2Points;

    

end