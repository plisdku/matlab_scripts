% Extrude  Representation of extruded polyhedron
%
% Usage:
%
% ex = Extrude('Vertices', xyVertFunc, 'Z', heightFunc, 'Permittivity',
%   epsString);
%
%   xyVertFunc:     function of parameter vector, returning x coordinates
%                   in the first column and y coordinates in the second
%                   column.  The xy pairs bound the region of material to
%                   fill in.
%
%   heightFunc:     function of parameter vector, returning a two-element
%                   column vector [z0 z1].  The xy polygon specified by the
%                   'Vertices' argument will be extruded vertically from
%                   z=z0 to z=z1.
%
%   Permittivity:   name of the permittivity model representing this
%                   material
%
%   Permeability:   name of the permeability model representing this
%                   material
%
classdef Extrude < t6.model.Node
    
    properties
        vertXFunc = @(p) [];
        vertYFunc = @(p) [];
        permittivity = 'none';
        permeability = 'none';
        zFunc = [];
    end
    
    methods
        
        % Extrude(points)
        %
        % points should be a function that provides an Nx2 array of x-y
        % coordinates of points.  Yeah!  All right!
        %
        % Extrude('Vertices', @(p) [0 0; 1 0; 0 1], 'Z', @(p) [0 1]', ...
        %   'Permittivity', 'Air', 'Permeability', 'Air');
        function obj = Extrude(varargin)
            if nargin > 0
                
                X.Vertices = [];
                X.Permittivity = '';
                X.Permeability = '';
                X.Z = [];
                X = parseargs(X, varargin{:});
                
                if ~isa(X.Vertices, 'function_handle')
                    error('Vertices must be a function handle');
                end
                if ~isa(X.Z, 'function_handle');
                    error('Z must be a function handle');
                end
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
                obj.zFunc = X.Z;
                
                firstCol = @(A) A(:,1);
                secondCol = @(A) A(:,2);
                obj.vertXFunc = @(p) firstCol(X.Vertices(p));
                obj.vertYFunc = @(p) secondCol(X.Vertices(p));
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
            xx = obj.vertXFunc(params);
            yy = obj.vertYFunc(params);
            zz = obj.zFunc(params);
            
            assert(size(xx,2) == 1);
            assert(size(yy,2) == 1);
            assert(size(xx,1) == size(yy,1));
            numVertices = size(xx,1);
            
            vertexTable = zeros(2*numVertices,1);
            vertexTable(:,1) = [xx; xx];
            vertexTable(:,2) = [yy; yy];
            vertexTable(1:numVertices,3) = zz(1);
            vertexTable(numVertices+1:2*numVertices,3) = zz(2);
            
            myVerts = zeros(6*numVertices, 1);
            myVerts(1:3:end) = vertexTable(:,1);
            myVerts(2:3:end) = vertexTable(:,2);
            myVerts(3:3:end) = vertexTable(:,3);
            
            % ---- Make a table of the actual jacobians (numbers too!)
            
            myJacobian = sparse(length(myVerts), size(params, 1));
            
            Dxx = jacobian(obj.vertXFunc, params);
            Dyy = jacobian(obj.vertYFunc, params);
            Dzz = jacobian(obj.zFunc, params);
            
            for vv = 1:numVertices
                row = 1 + 3*(vv-1);
                myJacobian(row,:) = Dxx(vv,:);
                myJacobian(row+1,:) = Dyy(vv,:);
                myJacobian(row+2,:) = Dzz(1,:);
            end
            
            for vv = numVertices+1:2*numVertices
                row = 1 + 3*(vv-1);
                myJacobian(row,:) = Dxx(vv-numVertices,:);
                myJacobian(row+1,:) = Dyy(vv-numVertices,:);
                myJacobian(row+2,:) = Dzz(2,:);
            end
            
            % ---- Work out the triangular faces.
            % We'll use a constrained Delaunay triangulation.
            
            profile = [xx yy];
            outerConstraint = [(1:numVertices)', [2:numVertices, 1]'];
            dt = DelaunayTri(profile, outerConstraint);
            inside = inOutStatus(dt);
            
            botTris = dt.Triangulation(inside,:);
            topTris = dt.Triangulation(inside,:) + numVertices;
            
            % The bottom triangles face the wrong way so permute their
            % vertices.
            botTris = botTris(:,[1 3 2]);
            
            % Now the side tris.
            sideTrisBot = [ (1:numVertices)', [2:numVertices 1]', ...
                numVertices + [2:numVertices 1]' ];
            sideTrisTop = [ (1:numVertices)', ...
                numVertices + [2:numVertices 1]', ...
                numVertices + (1:numVertices)' ];
            
            faces = [botTris; topTris; sideTrisBot; sideTrisTop];
            
            m = { Mesh(myVerts, faces, myJacobian, obj.permittivity,...
                obj.permeability) };
        end
        
    end
    
    
end