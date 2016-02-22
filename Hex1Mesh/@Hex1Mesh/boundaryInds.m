function [BoundVerts, BoundEle] = boundaryInds(T, varargin)
% Computing boundary elements and boundary indecies
% [BoundVerts, BoundEle]  = T.boundaryInds({dimension_1, value_1, dimension_2, value_2, ..., dimension_n, value_n })
% where dimension_i is e.g. 'x' or 'y' or 'z' and value_i is e.g. max(xnod) or min(znod)

if nargin < 2
    error('Not enough input')
end

if ~isempty(T.SurfaceP2Info)
    SurfEle = T.SurfaceP2Info.CutElements;
elseif ~isempty(T.SurfaceP1Info)
    SurfEle = T.SurfaceP1Info.CutElements;
else
    error('No surface exists!')
end

if ~iscell(varargin{1})
    error('input is not a cell!')
end

data = varargin{1};
dims = data(1:2:end);
vals = data(2:2:end);

nB = length(dims);

Bvert = [];
Bele = [];

xc = T.XC;
yc = T.YC;
zc = T.ZC;
nodes = T.Connectivity;



for i = 1:nB
    if length(vals{i}) > 1
        error('Content of value must be scalar')
    end
    
    if strcmpi(dims{i},'x')
        Bvert = [Bvert; find( xc >= vals{i}-eps*10 & xc <= vals{i}+eps*10 )];
        Bele = [Bele; find(sum(xc(nodes) >= vals{i}-eps*10 & xc(nodes) <= vals{i}+eps*10,2))];
    elseif strcmpi(dims{i},'y')
        Bvert = [Bvert; find( yc >= vals{i}-eps*10 & yc <= vals{i}+eps*10 )];
        Bele = [Bele; find(sum(yc(nodes) >= vals{1}-eps*10 & yc(nodes) <= vals{1}+eps*10,2))];
    elseif strcmpi(dims{i},'z')
        Bvert = [Bvert; find( zc >= vals{i}-eps*10 & zc <= vals{i}+eps*10 )];
        Bele = [Bele; find(sum(zc(nodes) >= vals{i}-eps*10 & zc(nodes) <= vals{i}+eps*10,2))];
    else
        error('dimension content must be ''x'' and/or ''y'' and/or ''z'' ')
    end
    
end



BoundEle = SurfEle(ismember(SurfEle,Bele));
Se = unique(nodes(BoundEle,:));
BoundVerts = Se(ismember(Se, Bvert ));

if find(strcmpi(varargin,'viz'), 1)
    T.vizMesh(BoundEle); hold on;
    plot3(xc(BoundVerts),yc(BoundVerts),zc(BoundVerts),'b*')
end

