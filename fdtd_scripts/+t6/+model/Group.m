classdef Group < t6.model.Node
    
    properties
        children = {};
    end
    
    methods
        function obj = Group(varargin)
            if nargin > 0
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
                    m{end+1} = childMeshes{mm};
                end
            end
        end
        
    end
end