classdef Rotate < t6.model.Node
    
    properties
        children = {};
        matrix = @(p) eye(3);
    end
    
    methods
        function obj = Rotate(matrix, varargin)
            if nargin > 0
                obj.matrix = matrix;
                obj.children = varargin;
            end
        end
        
        function m = meshes(obj, varargin)
            import t6.model.*
            
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            m = {};
            
            for cc = 1:length(obj.children)
                childMeshes = obj.children{cc}.meshes(params);
                for mm = 1:length(childMeshes)
                    v0 = childMeshes{mm}.vertices;
                    jac0 = childMeshes{mm}.jacobian;
                    
                    %numVerts = length(v0)/3;
                    
                    myVertices = blockProduct(...
                        obj.matrix(params), v0, [3 3], [3 1]);
                    
                    if ~isempty(params)
                        myJacobian = ...
                            blockProduct(obj.matrix(params), ...
                                jac0, [3 3], [3 1]) + ...
                            blockProduct(jacobian(obj.matrix, params), ...
                                v0, [3 3], [3 1]);
                    else
                        myJacobian = sparse(length(v0), 0);
                    end
                    
                    m{end+1} = Mesh(myVertices,...
                        childMeshes{mm}.faces,...
                        myJacobian, ...
                        childMeshes{mm}.permittivity,...
                        childMeshes{mm}.permeability);
                end
            end
        end
        
    end
    
end