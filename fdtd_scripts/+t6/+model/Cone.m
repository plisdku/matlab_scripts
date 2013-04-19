classdef Cone < t6.model.Node
% Cone   Representation of a cone!
%
% Constructor example:
%
% c = Cone('BaseX', @(p) [0 1 1 0], 'BaseY', @(p) [0 0 1 1], ...
%       'BaseZ', @(p) 0, 'TipX', @(p) 0.5, 'TipY', @(p) 0.5, ...
%       'TipZ', @(p) 1, 'Permittivity', 'Air', 'Permeability', 'Air');

    properties
        baseX = @(p) 0;
        baseY = @(p) 0;
        baseZ = @(p) 0;
        tipX = @(p) 0;
        tipY = @(p) 0;
        tipZ = @(p) 0;
        permittivity = 'none';
        permeability = 'none';
    end
    
    
    methods
        
        function obj = Cone(varargin)
            
            X.BaseX = [];
            X.BaseY = [];
            X.BaseZ = [];
            X.TipX = [];
            X.TipY = [];
            X.TipZ = [];
            X.Permittivity = '';
            X.Permeability = '';
            X = parseargs(X, varargin{:});
            
            obj.permittivity = X.Permittivity;
            obj.permeability = X.Permeability;
            obj.baseX = X.BaseX;
            obj.baseY = X.BaseY;
            obj.baseZ = X.BaseZ;
            obj.tipX = X.TipX;
            obj.tipY = X.TipY;
            obj.tipZ = X.TipZ;
        end
        
        function m = meshes(obj, varargin)
            import t6.model.*
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            row = @(A) reshape(A, 1, []);
            col = @(A) reshape(A, [], 1);
            bx = row(obj.baseX(params));
            by = row(obj.baseY(params));
            bz = row(obj.baseZ(params));
            tx = row(obj.tipX(params));
            ty = row(obj.tipY(params));
            tz = row(obj.tipZ(params));
            
            % The tip vertex is the last one!
            vertexTable = [ bx tx; by ty; bz(1)*ones(size(bx)) tz ];
            myVerts = vertexTable(:);
            
            m = vertexTable;
            
            % Jacobian table
            
            Dx = [jacobian(obj.baseX, params); ...
                jacobian(obj.tipX, params)];
            Dy = [jacobian(obj.baseY, params); ...
                jacobian(obj.tipY, params)];
            Dz = [repmat(jacobian(obj.baseZ, params), numel(bx), 1); ...
                jacobian(obj.tipZ, params)];
            
            myJacobian = sparse(3*size(Dx, 1), size(Dx,2));
            
            if numel(params) > 0
                myJacobian(1:3:end,:) = Dx;
                myJacobian(2:3:end,:) = Dy;
                myJacobian(3:3:end,:) = Dz;
            end
            
            % Faces!
            
            numSides = numel(bx);
            faces = zeros(numSides, 3);
            
            faces(:,1) = 1:numSides;
            faces(:,2) = [2:numSides, 1];
            faces(:,3) = numSides+1;
            
            % Add the bottom faces.
            
            profile = [bx' by'];
            outerConstraint = [(1:numSides)', [2:numSides,1]'];
            dt = DelaunayTri(profile, outerConstraint);
            inside = inOutStatus(dt);
            tris = dt.Triangulation(inside,:);
            
            faces = [faces; fliplr(tris)];
            
            m = { Mesh(myVerts, faces, myJacobian, obj.permittivity, ...
                obj.permeability) };
        end
    
    end
end


































