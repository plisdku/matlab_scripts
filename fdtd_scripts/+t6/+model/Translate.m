classdef Translate < t6.model.Node
% Translate(distance, children)
    
    properties
        children = {};
        dist = @(p) [0 0 0]';
    end
    
    methods
        function obj = Translate(distance, varargin)
            if nargin > 0
                obj.children = varargin;
                obj.dist = distance;
                
                assert(isa(obj.dist, 'function_handle'));
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
                    
                    numVerts = length(v0)/3;
                    
                    myDisplacement = repmat(obj.dist(params), numVerts, 1);
                    myJacobian = repmat(jacobian(obj.dist, params), ...
                        numVerts,1);
                    
                    m{end+1} = Mesh(v0 + reshape(myDisplacement',[],1),...
                        childMeshes{mm}.faces,...
                        jac0 + myJacobian, ...
                        childMeshes{mm}.permittivity,...
                        childMeshes{mm}.permeability);
                end
            end
        end
        
    end
end