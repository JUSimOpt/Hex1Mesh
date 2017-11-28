classdef HexP1MeshAbaqus < Hex1Mesh
    %HexP1MeshAbaqus
    %     8-----7
    %    /|    /|
    %   5-----6 |
    %   | 4...|.3
    %   |/    |/ 
    %   1-----2
    
    properties
        nodes
        P
        BC % Struct Boundary conditions %TODO: Remove
        NodeSets % contains Node sets
        ElementSets
        BoundarieConditions
        Loads        
    end
    
    methods
        function obj = HexP1MeshAbaqus(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
            obj = obj@Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
            map = [1,3,4,2,5,7,8,6];
            obj.Connectivity = obj.Connectivity(:,map);
            obj.nodes = obj.Connectivity;
            obj.P = obj.Points;
            
            
            
        end
        
        function [fi, B, detJ] = BasisFunctionsParametric(obj,ielement,xi,eta,zeta)
            iv = obj.nodes(ielement,:);
            XC = obj.P(iv,:);
            [fi, B, detJ] = baseHex1(XC, xi, eta, zeta);
            
        end
        function [fi, B,Bxi,Beta,Bzeta, detJ] = BasisFunctionsParametric1Point(obj,ielement,xi,eta,zeta)
            iv = obj.nodes(ielement,:);
            XC = obj.P(iv,:);
            [fi, B,Bxi,Beta,Bzeta, detJ] = baseHex1Point(XC);
        end
        
        function BasisFunctionsSpatial(ielement,x,y,z)
            error('Not implemented. TODO')
        end
        
    
    end
    
    methods(Static)
        function [GP,GW] = GaussQuadrature(order)
            [GP,GW] = Gauss3DHex(order);
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
        
        function baseHex1_inMidPoint(varargin)
            error('Not Implemented, use the other!!!')
        end
        function baseHex1_param(varargin)
            error('Not Implemented, use the other!!!')
        end
        function baseHexP1(varargin)
            error('Not Implemented, use the other!!!')
        end
        function baseHexP1_Parameter(varargin)
            error('Not Implemented, use the other!!!')
        end
        
    end
end

function [fi, B, detJ] = baseHex1(XC0, xi, eta, zeta)
    
% Unit cube with node one in 0,0,0
%     8-----7
%    /|    /|
%   5-----6 |
%   | 4...|.3
%   |/    |/ 
%   1-----2
%
% xi in [0,1]
% eta in [0,1]
% zeta in [0,1]


fi = [-(eta - 1)*(xi - 1)*(zeta - 1);...
        xi*(eta - 1)*(zeta - 1);...
             -eta*xi*(zeta - 1);...
        eta*(xi - 1)*(zeta - 1);...
        zeta*(eta - 1)*(xi - 1);...
             -xi*zeta*(eta - 1);...
                    eta*xi*zeta;...
             -eta*zeta*(xi - 1)].';
             

dfidxi = [-(eta - 1)*(zeta - 1);
  (eta - 1)*(zeta - 1);
       -eta*(zeta - 1);
        eta*(zeta - 1);
        zeta*(eta - 1);
       -zeta*(eta - 1);
              eta*zeta;
             -eta*zeta];
                  
dfideta = [-(xi - 1)*(zeta - 1);
        xi*(zeta - 1);
       -xi*(zeta - 1);
  (xi - 1)*(zeta - 1);
        zeta*(xi - 1);
             -xi*zeta;
              xi*zeta;
       -zeta*(xi - 1)];
                  
dfidzeta = [-(eta - 1)*(xi - 1);
        xi*(eta - 1);
             -eta*xi;
        eta*(xi - 1);
  (eta - 1)*(xi - 1);
       -xi*(eta - 1);
              eta*xi;
       -eta*(xi - 1)];


    X0 = XC0(:,1);
    Y0 = XC0(:,2);
    Z0 = XC0(:,3);

    J = [dfidxi.'*X0(:), dfidxi.'*Y0(:), dfidxi.'*Z0(:);...
         dfideta.'*X0(:), dfideta.'*Y0(:), dfideta.'*Z0(:);...
         dfidzeta.'*X0(:), dfidzeta.'*Y0(:), dfidzeta.'*Z0(:)];

    detJ = det(J); 

    B = J\[dfidxi.';dfideta.';dfidzeta.'];

end

function [fi, B, Bxi, Beta, Bzeta, detJ] = baseHex1Point(XC0)
    
% Unit cube with node one in 0,0,0
%     8-----7
%    /|    /|
%   5-----6 |
%   | 4...|.3
%   |/    |/ 
%   1-----2
%
% xi in [0,1]
% eta in [0,1]
% zeta in [0,1]


% fi = [-(eta - 1)*(xi - 1)*(zeta - 1);...
%         xi*(eta - 1)*(zeta - 1);...
%              -eta*xi*(zeta - 1);...
%         eta*(xi - 1)*(zeta - 1);...
%         zeta*(eta - 1)*(xi - 1);...
%              -xi*zeta*(eta - 1);...
%                     eta*xi*zeta;...
%              -eta*zeta*(xi - 1)].';            

% dfidxi = [-(eta - 1)*(zeta - 1);
%   (eta - 1)*(zeta - 1);
%        -eta*(zeta - 1);
%         eta*(zeta - 1);
%         zeta*(eta - 1);
%        -zeta*(eta - 1);
%               eta*zeta;
%              -eta*zeta].';
%                   
% dfideta = [-(xi - 1)*(zeta - 1);
%         xi*(zeta - 1);
%        -xi*(zeta - 1);
%   (xi - 1)*(zeta - 1);
%         zeta*(xi - 1);
%              -xi*zeta;
%               xi*zeta;
%        -zeta*(xi - 1)].';
%                   
% dfidzeta = [-(eta - 1)*(xi - 1);
%         xi*(eta - 1);
%              -eta*xi;
%         eta*(xi - 1);
%   (eta - 1)*(xi - 1);
%        -xi*(eta - 1);
%               eta*xi;
%        -eta*(xi - 1)].';

    fi = [0.1250    0.1250    0.1250    0.1250    0.1250    0.1250    0.1250    0.1250]; 
    dfidxi = [-0.2500    0.2500    0.2500   -0.2500   -0.2500    0.2500    0.2500   -0.2500];
    dfideta = [-0.2500   -0.2500    0.2500    0.2500   -0.2500   -0.2500    0.2500    0.2500];
    dfidzeta = [-0.2500   -0.2500   -0.2500   -0.2500    0.2500    0.2500    0.2500    0.2500];


    X0 = XC0(:,1);
    Y0 = XC0(:,2);
    Z0 = XC0(:,3);

    J = [dfidxi*X0(:), dfidxi*Y0(:), dfidxi*Z0(:);...
         dfideta*X0(:), dfideta*Y0(:), dfideta*Z0(:);...
         dfidzeta*X0(:), dfidzeta*Y0(:), dfidzeta*Z0(:)];

    detJ = det(J); 

%     B = J\[dfidxi;dfideta;dfidzeta];
    Bh = [dfidxi;dfideta;dfidzeta];

    Bhxi = [0,    0,    0,    0,    0,    0,   0,    0;
            0.5, -0.5,  0.5, -0.5,  0.5, -0.5, 0.5, -0.5;
            0.5, -0.5, -0.5,  0.5, -0.5,  0.5, 0.5, -0.5];
        
    Bheta = [0.5, -0.5,  0.5, -0.5,  0.5, -0.5, 0.5, -0.5;
             0,    0,    0,    0,    0,    0,   0,    0;
             0.5,  0.5, -0.5, -0.5, -0.5, -0.5, 0.5,  0.5];
          
    Bhzeta = [0.5, -0.5, -0.5,  0.5, -0.5,  0.5, 0.5, -0.5;
              0.5,  0.5, -0.5, -0.5, -0.5, -0.5, 0.5,  0.5;
              0,    0,    0,    0,    0,    0,   0,    0];
        
    BB = J\[Bh,Bhxi,Bheta,Bhzeta];
    B = BB(:,1:8);
    Bxi = BB(:,9:16);
    Beta = BB(:,17:24);
    Bzeta = BB(:,25:32);
    
    
end

function [GP,GW] = Gauss3DHex(order)
    [x,w] = gauss(order, 0,1);
    %x in [0,1]
    [X,Y,Z]=meshgrid(x,x,x);
    GP = [X(:),Y(:),Z(:)];
    
    [WX,WY,WZ]=meshgrid(w,w,w);
    GW = [WX(:),WY(:),WZ(:)];
    
    GW = prod(GW,2);
end


function [x, w, A] = gauss(n, a, b)

%------------------------------------------------------------------------------
% gauss.m
%------------------------------------------------------------------------------
%
% Purpose:
%
% Generates abscissas and weigths on I = [ a, b ] (for Gaussian quadrature).
%
%
% Syntax:
%
% [x, w, A] = gauss(n, a, b);
%
%
% Input:
%
% n    integer    Number of quadrature points.
% a    real       Left endpoint of interval.
% b    real       Right endpoint of interval.
%
%
% Output:
%
% x    real       Quadrature points.
% w    real       Weigths.
% A    real       Area of domain.
%------------------------------------------------------------------------------


% 3-term recurrence coefficients:
n = 1:(n - 1);
beta = 1 ./ sqrt(4 - 1 ./ (n .* n));

% Jacobi matrix:
J = diag(beta, 1) + diag(beta, -1); 


% Eigenvalue decomposition:

%
% e-values are used for abscissas, whereas e-vectors determine weights.
%

[V, D] = eig(J);
x = diag(D);


% Size of domain:
A = b - a;


% Change of interval:

%
% Initally we have I0 = [ -1, 1 ]; now one gets I = [ a, b ] instead.
%
% The abscissas are Legendre points.
%

if ~(a == -1 && b == 1)
  x = 0.5 * A * x + 0.5 * (b + a);
end


% Weigths:
w = V(1, :) .* V(1, :);
end

