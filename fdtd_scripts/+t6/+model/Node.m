classdef Node < handle
    
    methods (Abstract)
        m = meshes(obj, varargin); % vertices, faces, jacobians, material
    end
    
end