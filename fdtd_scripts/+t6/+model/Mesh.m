classdef Mesh < handle
    
    properties
        vertices
        faces
        permittivity
        permeability
        jacobian
    end
    
    methods
        
        function obj = Mesh(v, f, j, permittivity, permeability)
            if nargin > 0 % default constructor has same signature, sigh
                obj.vertices = v;
                obj.faces = f;
                obj.jacobian = j;
                obj.permittivity = permittivity;
                obj.permeability = permeability;
            end
        end
        
        function f = freeDirections(obj, varargin)
            f = zeros(size(obj.vertices,1)/3, 3);
            
            for vv = 1:size(f,1)
                for xyz = 1:3
                    if any(obj.jacobian(3*(vv-1)+xyz,:))
                        f(vv,xyz) = 1;
                    end
                end
            end
            
        end
        
        function v = patchVertices(obj, varargin)
            vertexColumn = obj.vertices(varargin{:});
            v = [vertexColumn(1:3:end), vertexColumn(2:3:end),...
                vertexColumn(3:3:end)];
        end
        
        function b = bounds(obj, varargin)
            
            v = obj.patchVertices(varargin{:});
            minCoord = min(v);
            maxCoord = max(v);
            
            b = [minCoord maxCoord];
        end
    end
    
end