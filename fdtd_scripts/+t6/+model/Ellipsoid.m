% Ellipsoid    Representation of an ellipsoid
%
% Constructor example: create a sphere
%
% s = Ellipsoid(@(p) [0 0 0 1 1 1], 'NumThetas', 10, 'NumPhis', 10, ...
%   'Permittivity', 'Air', 'Permeability', 'Air');
%

classdef Ellipsoid < t6.model.Node
    
    properties
        numThetas = 10;
        numPhis = 10;
        bounds = @(params) [0 0 0 1 1 1]';
        centerFunc = @(params) [0 0 0]';
        radiusFunc = @(params) [1 1 1]';
        permittivity = 'none';
        permeability = 'none';
    end
    
    methods
    
        function obj = Ellipsoid(bounds, varargin)
            if nargin > 0
                X.NumThetas = 10;
                X.NumPhis = 10;
                X.Permittivity = '';
                X.Permeability = '';
                X = parseargs(X, varargin{:});
                
                if X.NumThetas < 3
                    error('NumThetas must be at least 3');
                end
                
                if X.NumPhis < 3
                    error('NumPhis must be at least 3');
                end
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
                obj.numThetas = X.NumThetas;
                obj.numPhis = X.NumPhis;
                obj.bounds = @(p) transpose(bounds(p));
                
                rectCenter = @(r) 0.5*(r(4:6)+r(1:3));
                rectRadius = @(r) 0.5*(r(4:6)-r(1:3));
                obj.centerFunc = @(p) rectCenter(obj.bounds(p));
                obj.radiusFunc = @(p) rectRadius(obj.bounds(p));
            end
        end
            
        function m = meshes(obj, varargin)
            import t6.model.*
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            thetas = linspace(0, pi, obj.numThetas);
            phis = linspace(0, 2*pi, obj.numPhis+1);
            phis = phis(1:end-1);
            
            % ---- Make a table of the vertex positions (numbers not funcs!)
            r = obj.radiusFunc(params);
            c = obj.centerFunc(params);
            
            vertexTable = zeros(2 + obj.numPhis * (obj.numThetas-2));
            
            % Index functions
            vTop = 1;
            vBot = 2;
            vSide = @(th, ph) vBot + mod(ph-1, obj.numPhis) + 1 + ...
                (th-2)*obj.numPhis;
            
            numVerts = 2 + obj.numPhis*(obj.numThetas-2);
            vThetas = zeros(numVerts,1);
            vPhis = zeros(numVerts,1);
            
            vThetas([vTop vBot]) = [0 pi];
            vPhis([vTop vBot]) = 0;
            
            for pp = 1:obj.numPhis
            for tt = 2:obj.numThetas-1
                vThetas(vSide(tt,pp)) = thetas(tt);
                vPhis(vSide(tt,pp)) = phis(pp);
            end
            end
            
            % Vertex positions
            myVerts = zeros(3*numVerts,1);
            myVerts(1:3:end) = c(1) + r(1)*sin(vThetas).*cos(vPhis);
            myVerts(2:3:end) = c(2) + r(2)*sin(vThetas).*sin(vPhis);
            myVerts(3:3:end) = c(3) + r(3)*cos(vThetas);
            
            % ---- Make a table of the vertex positions (numbers not funcs!)
            
            myJacobian = sparse(length(myVerts), size(params, 1));
            
            Dr = jacobian(obj.radiusFunc, params);
            Dc = jacobian(obj.centerFunc, params);
            
            for vv = 1:numVerts
                row = 1 + 3*(vv-1);
                theta = vThetas(vv);
                phi = vPhis(vv);
                myJacobian(row,:) = Dc(1,:) + Dr(1,:)*sin(theta)*cos(phi);
                myJacobian(row+1,:) = Dc(2,:) + Dr(2,:)*sin(theta)*sin(phi);
                myJacobian(row+2,:) = Dc(3,:) + Dc(3,:)*cos(theta);
            end
            
            %% ---- Work out the triangular faces
            numTopTris = obj.numPhis;
            numBotTris = numTopTris;
            numSideTris = (obj.numThetas-2)*obj.numPhis;
            numFaces = numTopTris + numBotTris + numSideTris;
            
            faces = zeros(numFaces, 3);
            
            % Top tris
            topFaces = [];
            topFaces(1:numTopTris, 1) = vTop;
            topFaces(1:numTopTris, 2) = vSide(2, 1:obj.numPhis);
            topFaces(1:numTopTris, 3) = vSide(2, 1 + (1:obj.numPhis));
            
            % Side tris
            sideFaces = [];
            for pp = 1:obj.numPhis
            for tt = 2:obj.numThetas-2
                sideFaces(end+1,:) = [vSide(tt,pp+1) vSide(tt,pp) vSide(tt+1,pp)];
                sideFaces(end+1,:) = [vSide(tt,pp+1) vSide(tt+1,pp) vSide(tt+1,pp+1)];
            end
            end
            
            % Bottom tris
            bottomFaces = [];
            bottomFaces(1:numBotTris, 1) = vBot;
            bottomFaces(1:numBotTris, 2) = vSide(obj.numThetas-1, 1:obj.numPhis);
            bottomFaces(1:numBotTris, 3) = vSide(obj.numThetas-1, 1 + (1:obj.numPhis));
            
            faces = [topFaces; sideFaces; bottomFaces];
            
            m = { Mesh(myVerts, faces, myJacobian, obj.permittivity, obj.permeability) };
        end
    end % methods
end


