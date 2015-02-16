function hp = vizP1Surf(T, varargin)
%VIZP1SURF Visualize P1 surface
%   hp = vizP1Surf(T)
%   hp = vizP1Surf(T, parameters)
%   parameters:
%   showMesh

if isempty(T.SurfaceP1)
    error('Surface must exist!')
end

if ~isempty(T.SurfaceP1Triangulation)
    if nargin > 1
    	hp.fig = xfigure(varargin{1});
    else
        hp.fig = xfigure;
    end

    axis equal; hold on;
    FV.Vertices = T.SurfaceP1Points;
    FV.Faces = T.SurfaceP1Triangulation;
    shading interp
    hp.light = light;
    hp.patch = patch(FV,'FaceColor','c','FaceLighting','gouraud');
    [az,el]=view(136,16);
    hp.az = az;
    hp.el = el;
    
    param = 'showMesh';
    if isenabled(param,varargin)
        hp.mesh = showMesh(T);
    end
    
    

else
    surfX = T.SurfaceP1Points;
    xx = reshape(surfX(:,1),3,[]);
    yy = reshape(surfX(:,2),3,[]);
    zz = reshape(surfX(:,3),3,[]);
    hp.fig = xfigure; 
    if nargin > 1
    	hp.fig = xfigure(varargin{1});
    else
        hp.fig = xfigure;
    end
    
    axis equal; hold on;
    hp.patch = patch(xx,yy,zz,'r');
    [az,el]=view(130,14);
    hp.az = az;
    hp.el = el;
    
    param = 'showMesh';
    if isenabled(param,varargin)
        hp.mesh = showMesh(T);
    end
end

ntri = length(T.SurfaceP1Points)/3;
nCutEle = T.SurfaceP1Info.NCutElements;

title([num2str(ntri),' triangles on ',num2str(nCutEle),' cut Hex1 elements'])

end





function mesh = showMesh(T)
    ele = 1:size(T.Connectivity,1);
    ele = ele(:);
    fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
    mesh = patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
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