classdef Rect < t6.model.Node
% Rect  Representation of an axis-aligned rectangle
% 
% Constructor example:
%
% r = Rect(@(p) [0 0 0 1 1 1], 'Permittivity', 'Air');
%
% The bounding box should be a column vector.  (WHY IS THIS?)
% 
    properties
        func = @(params) [0 0 0 1 1 1]';
        permittivity = 'none';
        permeability = 'none';
    end
    
    methods
        
        function obj = Rect(func, varargin)
            if nargin > 0
                obj.func = func;
                
                X.Permittivity = '';
                X.Permeability = '';
                X = parseargs(X, varargin{:});
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
            end
        end
        
        function m = meshes(obj, varargin)
            import t6.model.*
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            faces = [...
                1 6 5; ...
                1 2 6; ...
                4 7 8; ...
                4 3 7; ...
                1 4 2; ...
                1 3 4; ...
                5 8 7; ...
                5 6 8; ...
                1 7 3; ...
                1 5 7; ...
                2 8 6; ...
                2 4 8; ...
            ];
            
            corners = transpose(obj.func(params));
            Dcorners = t6.model.jacobian(@(p) transpose(obj.func(p)), params);
            
            myVerts = zeros(24,1);
            myJacobian = sparse(24, size(params,1));
            
            row = 1;
            for zz = 0:1, for yy = 0:1, for xx = 0:1
                myVerts(row:row+2,:) = corners([1+3*xx, 2+3*yy, 3+3*zz]);
                
                myJacobian(row:row+2,:) = Dcorners(...
                    [1+3*xx, 2+3*yy, 3+3*zz], :);
                
                row = row + 3;
            end,end,end
            
            m = { Mesh(myVerts, faces, myJacobian, obj.permittivity,...
                obj.permeability) };
            
        end
        
    end
    
    
end