classdef Cylinder < t6.model.Node
    
    properties
        numSides = 3;
        bounds = @(params) [0 0 0 1 1 1]';
        centerFunc = @(params) [0 0 0]';
        radiusFunc = @(params) [1 1 1]';
        permittivity = 'none';
        permeability = 'none';
    end
    
    methods
        
        function obj = Cylinder(bounds, varargin)
            if nargin > 0
                
                X.NumSides = 3;
                X.Permittivity = '';
                X.Permeability = '';
                X = parseargs(X, varargin{:});
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
                obj.numSides = X.NumSides;
                obj.bounds = bounds;
                
                rectCenter = @(r) 0.5*(r(4:6) + r(1:3));
                rectRadius = @(r) 0.5*(r(4:6) - r(1:3));
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
            
            % ---- Make a table of the actual vertex positions (numbers!)
            r = obj.radiusFunc(params);
            c = obj.centerFunc(params);
            boundingBox = obj.bounds(params);
            
            vertexTable = zeros(2*obj.numSides,1);
            
            for vv = 1:obj.numSides
                theta = 2*pi*vv/obj.numSides;
                vertexTable(vv, 1) = c(1) + r(1)*cos(theta);
                vertexTable(vv, 2) = c(2) + r(2)*sin(theta);
                vertexTable(vv, 3) = boundingBox(3);
            end
            
            for vv = 1:obj.numSides
                theta = 2*pi*vv/obj.numSides;
                vertexTable(vv+obj.numSides, 1) = c(1) + r(1)*cos(theta);
                vertexTable(vv+obj.numSides, 2) = c(2) + r(2)*sin(theta);
                vertexTable(vv+obj.numSides, 3) = boundingBox(6);
            end
            
            myVerts = zeros(6*obj.numSides, 1);
            myVerts(1:3:end) = vertexTable(:,1);
            myVerts(2:3:end) = vertexTable(:,2);
            myVerts(3:3:end) = vertexTable(:,3);
            
            % ---- Make a table of the actual jacobians (numbers too!)
            
            myJacobian = sparse(length(myVerts), size(params, 1));
            
            % Derivative of vertex position:
            % vx: Dc(1) + Dr(1)*cos(theta)
            % vy: Dc(2) + Dr(2)*sin(theta)
            % vz: DboundingBox(3) or DboundingBox(6)
            
            Dr = jacobian(obj.radiusFunc, params);
            Dc = jacobian(obj.centerFunc, params);
            Dbb = jacobian(obj.bounds, params);
            
            for vv = 1:obj.numSides
                row = 1 + 3*(vv-1);
                theta = 2*pi*vv/obj.numSides;
                myJacobian(row,:) = Dc(1,:) + Dr(1,:)*cos(theta);
                myJacobian(row+1,:) = Dc(2,:) + Dr(2,:)*sin(theta);
                myJacobian(row+2,:) = Dbb(3,:);
            end
            
            for vv = obj.numSides+1:2*obj.numSides
                row = 1 + 3*(vv-1);
                theta = 2*pi*vv/obj.numSides;
                myJacobian(row,:) = Dc(1,:) + Dr(1,:)*cos(theta);
                myJacobian(row+1,:) = Dc(2,:) + Dr(2,:)*sin(theta);
                myJacobian(row+2,:) = Dbb(6,:);
            end
            
            % ---- Work out the triangular faces.
            faces = zeros(4*obj.numSides - 4, 3);
            
            % sides, bottom edge first
            firstFace = 1;
            lastFace = obj.numSides;
            
            faces(firstFace:lastFace, 1) = 1:obj.numSides;
            faces(firstFace:lastFace, 2) = [2:obj.numSides, 1];
            faces(firstFace:lastFace, 3) = (1:obj.numSides)+obj.numSides;
            
            % sides, top edges
            firstFace = obj.numSides + 1;
            lastFace = 2*obj.numSides;
            
            faces(firstFace:lastFace, 1) = (1:obj.numSides)+obj.numSides;
            faces(firstFace:lastFace, 2) = [2:obj.numSides, 1];
            faces(firstFace:lastFace, 3) = [2:obj.numSides, 1] + obj.numSides;
            
            % bottom faces
            firstFace = 2*obj.numSides + 1;
            lastFace = firstFace + obj.numSides - 3;
            
            faces(firstFace:lastFace, 1) = 1;
            faces(firstFace:lastFace, 2) = 3:obj.numSides;
            faces(firstFace:lastFace, 3) = 2:(obj.numSides-1);
            
            % top face
            firstFace = lastFace + 1;
            lastFace = firstFace + obj.numSides - 3;
            
            faces(firstFace:lastFace, 1) = obj.numSides + 1;
            faces(firstFace:lastFace, 2) = (2:(obj.numSides-1)) + obj.numSides;
            faces(firstFace:lastFace, 3) = (3:obj.numSides) + obj.numSides;
            
            m = { Mesh(myVerts, faces, myJacobian, obj.permittivity,...
                obj.permeability) };
        end
        
    end
    
    
end