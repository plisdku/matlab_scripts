classdef AssignParameters < t6.model.Node
% AssignParameters(paramIndices, children)
    
    properties
        children = {}
        paramIndices = [];
    end
    
    methods
        function obj = AssignParameters(inIndices, varargin)
            if nargin > 0
                obj.children = varargin;
                obj.paramIndices = inIndices;
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
                childMeshes = obj.children{cc}.meshes(params(obj.paramIndices));
                
                for mm = 1:length(childMeshes)
                    
                    reducedJacobian = childMeshes{mm}.jacobian;
                    
                    fullJacobian = sparse(size(reducedJacobian, 1), ...
                        numel(params));
                    
                    fullJacobian(:, obj.paramIndices) = reducedJacobian;
                    
                    m{end+1} = Mesh(childMeshes{mm}.vertices, ...
                        childMeshes{mm}.faces, ...
                        fullJacobian, ...
                        childMeshes{mm}.permittivity, ...
                        childMeshes{mm}.permeability);
                end
            end
        end
    end
end
