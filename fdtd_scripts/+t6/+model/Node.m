classdef Node < handle
    
    methods (Abstract)
        m = meshes(obj, varargin); % vertices, faces, jacobians, material
    end
    
    methods
        function jj = jacobian(obj, parameters)
            jj = [];
            
            myMeshes = obj.meshes(parameters);
            
            for mm = 1:length(myMeshes)
                jj = [jj; myMeshes{mm}.jacobian];
            end
        end
        
        function vv = vertices(obj, parameters)
            vv = [];
            
            myMeshes = obj.meshes(parameters);
            
            for mm = 1:length(myMeshes)
                vv = [vv; myMeshes{mm}.vertices];
            end
        end
    end
end