classdef FreeMesh < t6.model.Node
    
    properties
        vertices = [];
        faces = [];
        freeDirections = [];
        permittivity = 'none';
        permeability = 'none';
        matrix = [];
    end
    
    
    methods
        
        function obj = FreeMesh(varargin)
            if nargin > 0
                X.Vertices = [];
                X.FreeDirections = [];
                X.Matrix = [];
                X.Faces = [];
                X.Permittivity = [];
                X.Permeability = [];
                X = parseargs(X, varargin{:});
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
                obj.vertices = X.Vertices;
                obj.matrix = X.Matrix;
                
                if ~isempty(X.Matrix)
                    X.FreeDirections = reshape(any(X.Matrix,2), 3, [])';
                end
                
                if ~isempty(X.FreeDirections)
                    obj.freeDirections = logical(X.FreeDirections);
                else
                    obj.freeDirections = false(size(obj.vertices));
                end
                
                obj.faces = X.Faces;
                
                assert(isequal(size(obj.freeDirections), ...
                    size(obj.vertices)));
            end
        end
        
        function m = meshes(obj, varargin)
            import t6.model.*
            if nargin > 1
                params = varargin{1};
            else
                params = zeros(nnz(obj.freeDirections), 1);
            end
            
            if ~isempty(obj.matrix)
                [myVerts, ~, myJacobian] = obj.matrixMesh(params);
            else
                [myVerts, ~, myJacobian] = obj.defaultMesh(params);
            end
            
            m = { Mesh(myVerts, obj.faces, myJacobian, obj.permittivity, ...
                obj.permeability) };
            
        end
        
        function [myVerts, myFreeDirs, myJacobian] = defaultMesh(obj, params)
            
            numParams = numel(params);
            if nnz(obj.freeDirections)
                assert(numParams == nnz(obj.freeDirections));
            end
            
            
            
            % All x first, then all y, then all z.
            movedVerts = obj.vertices;
            movedVerts(obj.freeDirections) = ...
                movedVerts(obj.freeDirections) + params;
            
            numVertices = size(obj.vertices, 1);
            myVerts = zeros(3*numVertices,1);
            myVerts(1:3:end) = movedVerts(:,1);
            myVerts(2:3:end) = movedVerts(:,2);
            myVerts(3:3:end) = movedVerts(:,3);
            
            myFreeDirs = zeros(3*numVertices,1);
            myFreeDirs([1:3:end, 2:3:end, 3:3:end],1) = ...
                obj.freeDirections(:);
            
            % Internally, vertices is a colvec [x0 y0 z0 x1 y1 z1]...
            % Parameters go to all x, all y,... [1 ?  ?  2  ?  ?  ] ...
            % I need to extract the numbers 1, 2 etc. and also their index
            % into the vertex colvec.  The index I'll call "coord" below.
            whichParams = 0*obj.freeDirections;
            whichParams(obj.freeDirections) = 1:numParams;
            whichParams = transpose(whichParams);
            [coord,~,param] = find(whichParams(:));
            
            % Rows of Jacobian are "coord" values, x0 y0 z0 x1 y1 z1 ...
            % Columns of Jacobian are parameter numbers ("param" above).
            % Values are all 1 for this FreeMesh.
            myJacobian = sparse(coord, param, ones(numParams,1), ...
                3*numVertices, numParams);
            
        end
        
        function [myVerts, myFreeDirs, myJacobian] = matrixMesh(obj, params)
            
            v0 = transpose(obj.vertices);
            myVerts = v0(:) + obj.matrix * params;
            myFreeDirs = any(obj.matrix, 2);
            myJacobian = obj.matrix;
            
        end
        
    end
    
    
end