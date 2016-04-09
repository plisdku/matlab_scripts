classdef Mesh < handle
    % Mesh object
    
    properties
        vertices
        faces
        permittivity
        permeability
        jacobian
    end
    
    methods
        
        function obj = Mesh(v, f, j, permittivity, permeability)
            % mesh = Mesh(vertices, faces, Jacobian, epsilon, mu)
            % 
            if nargin >= 2 % default constructor has same signature, sigh
                
                if size(v,2) == 3
                    v = v';
                    obj.vertices = v(:);
                elseif size(v,2) == 1              
                    obj.vertices = v;
                else
                    error('Vertex array must be Mx1 or Nx3');
                end
                obj.faces = f;
            end
            
            if nargin >= 3
                obj.jacobian = j;
            end
            
            if nargin >= 4
                obj.permittivity = permittivity;
            end
            
            if nargin >= 5
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