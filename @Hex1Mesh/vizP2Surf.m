function hv = vizP2Surf(T,varargin)
%VIZP2SURF Visualize P2 surface
%   hv = vizP2Surf(T, parameters)
%   Requires PlotP2
%   parameters:
%   nP2Ele: Number of subtriangles per P2 triangle
%   EdgeColor
%   FaceColor

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
    
    T.SurfaceP2Triangulation = tri;
    T.SurfaceP2Points = X;
end

tri = T.SurfaceP2Triangulation ;
X = T.SurfaceP2Points;

%% Default parameters
EdgeColor = 'none';
FaceColor = 'c';
nP2Ele = 4;

%% Optional parameters
param = 'EdgeColor';
if isenabled(param,varargin)
    EdgeColor = getoption(param,varargin);
end

param = 'FaceColor';
if isenabled(param,varargin)
    FaceColor = getoption(param,varargin);
end

param = 'nP2Ele';
if isenabled(param,varargin)
    nP2Ele = getoption(param,varargin);
end

%% P2 surface
if nargin > 1
    hv = PlotP2(tri,X,nP2Ele,'m',varargin{1});
else
    hv = PlotP2(tri,X,nP2Ele,'m');
end
hv.EdgeColor = EdgeColor;
hv.FaceColor = FaceColor;
view(149,31)

end

function rv = isenabled(mode, varargin)
%   ISENABLED  Checks if mode exists in the cell-array varargin.
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

function param = getoption(mode, varargin)
%   GETOPTION  Gets the parameter for the option 'mode'
%
%   getoption(mode,varargin) returns the parameter.
%   example:
%
%          varargin = {'contour', 'grid', 'view', [20,20]};
%          getoption('view',varargin)
%          ans =
%               20   20
%
%   NOTES: Input verification/ errorahandling should be added in the caller function to
%   keep things stable.
%
varargin = varargin{:};
ind1 = find(strcmpi(varargin,mode), 1);
if ~isempty(ind1)
    %Errorhandling    
    if ~iscell(varargin)
        varargin = {varargin};
    end
    if ind1+1 <= length(varargin)
        param = varargin{ind1+1};
    else
%         error(['No options are followed by the property ''', mode,''' '])
        param = [];
    end
else
    param = [];
end
end
