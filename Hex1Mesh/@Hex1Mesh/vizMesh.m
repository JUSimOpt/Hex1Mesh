function h = vizMesh(T,varargin)
%vizMesh Visualize the mesh in a figure
%
%   h = vizMesh()
%   h = vizMesh(figNumber)
%   h = vizMesh(figNumber,ele,properties)
%   ele is a list of elements to display
%   h contains visualization properties
%   figNumber sets the figure number
%   
%   Requirements: xfigure
%
%   properties:
%   'NodeNumbers'       -   Display node numbers
%   'ElementNumbers'    -   Display element numbers
%
%   Usage:
%   T.vizMesh( ... )
%   or
%   vizMesh(T, ... )
%   T is the Hex1Mesh class


fn = NaN;
if nargin > 1
    if isa(varargin{1},'double')
        fn = varargin{2};
    end
end
ele = 1:size(T.Connectivity,1);
if nargin > 2
    if isa(varargin{1},'double')
        ele = varargin{1};
    end
end
if ~ishold 
    if isnan(fn)
        if exist('xfigure','file') == 2
            h.fig = xfigure;
        else
            RequiredFileMissing('xfigure', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure.m')
            RequiredFileMissing('xfigure_KPF', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure_KPF.m')
            h.fig = xfigure;
        end
    else
        if exist('xfigure','file') == 2
            h.fig = xfigure(fn);
        else
            RequiredFileMissing('xfigure', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure.m')
            RequiredFileMissing('xfigure_KPF', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure_KPF.m')
            h.fig = xfigure(fn);
        end
    end
else
    h.fig = gcf;
end

ele = ele(:);
fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];

h.patch = patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
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
        h.NodeText = [h.NodeText; text(T.XC(i),T.YC(i),T.ZC(i),num2str(i),'BackgroundColor','w') ];
    end
end

if isenabled('ElementNumbers',varargin)
    if T.nele > 100
        warning('Cannot draw ElementNumbers, too many elements')
        return
    end
    
    h.EleText = [];
    for i = sort(ele)'
        xm = mean(T.XC(T.Connectivity(i,:)));
        ym = mean(T.YC(T.Connectivity(i,:)));
        zm = mean(T.ZC(T.Connectivity(i,:)));
        h.EleText = [h.EleText; text(xm,ym,zm,num2str(i),'BackgroundColor','y')];
    end
    
end


end

function RequiredFileMissing(filename, RemoteDestination)
    %If we're going to download a whole bunch of files, it is better to set
    % RequiredFilesDir to be a global and not have to ask the user to
    % specify a destination folder for every file...
    global RequiredFilesDir
    
    
    disp([filename,' is missing!'])
    disp(['Trying to download ',filename,' from ',RemoteDestination])
    
    
    if isempty(RequiredFilesDir)
        scriptPath = mfilename('class');
        [ScriptDir,~,~] = fileparts(scriptPath);
        DestDir = uigetdir(ScriptDir,['Select where to save ',filename,'. Make sure its either the script directory or a directory on the Path.']);
        if DestDir == 0
            error(['Failed to select folder, failed to install ',filename])
        end
        
        RequiredFilesDir = DestDir;
    end
    DestFile= [RequiredFilesDir,'/',filename];
    
    % Download the RemoteFile and save it to DestFile
    websave(DestFile,RemoteDestination);
    
    % Give up to 10 seconds for the file to show up, otherwise send error
    % message.
    tic
    while 1
        if exist('xfigure','file') == 2
            break
        end
        pause(0.1)
        t1 = toc;
        if t1 > 10
            error(['Failed to download ',filename,'! Timeout.'])
        end
    end

    
end


function rv = isenabled(mode, varargin)
    %   ISENABLED  Checks if mode exists in the cell-array varargin.
    %
    %   isenabled(mode,varargin{:}) return true or false.
    %   example:
    %
    %          varargin = {'Viz', 'ElementNumber', 'debug', [20,20]};
    %          isenabled('debug',varargin)
    %          ans =
    %               1
    %
    %   Author: Mirza Cenanovic (mirza.cenanovic@jth.hj.se)
    %   Date: 2013-05-02
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




