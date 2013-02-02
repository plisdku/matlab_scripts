% Heightmap   Representation of a height map, z = z(x,y)
%
% Usage:
%   h = Heightmap('X', @(p) linspace(0,1), 'Y', @(p) linspace(0,1), ...
%       'Z', @(p) peaks(100), 'Zbot', @(p) -1, ...
%       'Permittivity', 'Air', 'Permeability', 'Air');
% 
% Named parameters:
%       X   Function returning ascending array of x coordinates
%       Y   Function returning ascending array of y coordinates
%       Z   Function returning height of surface over each point, indexed (x,y)
%       Zbot Function returning the z coordinate of the bottom of the heightmap
%       Permittivity, Permeability  Names of material models to apply
%
% Heightmaps are always oriented such that "up" is the +z direction.  Use a
% Rotate object to re-orient them to other directions.
%
% If the top surface of the heightmap ever touches the bottom surface (i.e. if
% any(z(p) < zbot(p))) then the triangulator will crash.  To make a heightmap
% with holes in the bottom you will need to chop off the bottom with another 
% block of material (e.g. a Rect).
%
classdef Heightmap < t6.model.Node
    
    
    properties 
        permittivity = '';
        permeability = '';
        xFunc = @(p) 1;
        yFunc = @(p) 1;
        zFunc = @(p) 1;
        zBotFunc = @(p) 0;
    end
    
    
    methods
        function obj = Heightmap(varargin)
            
            if nargin > 0
                
                X.X = [];
                X.Y = [];
                X.Z = [];
                X.ZBot = [];
                X.Permittivity = '';
                X.Permeability = '';
                X = parseargs(X, varargin{:});
                
                obj.permittivity = X.Permittivity;
                obj.permeability = X.Permeability;
                obj.xFunc = X.X;
                obj.yFunc = X.Y;
                obj.zFunc = X.Z;
                obj.zBotFunc = X.ZBot;
                
            end
            
        end
        
        function [v nx ny] = vertices(obj, params)
            
            x = obj.xFunc(params);
            y = obj.yFunc(params);
            zz = obj.zFunc(params);
            
            nx = numel(x);
            ny = numel(y);
            
            [xx yy] = ndgrid(x, y);
            
            colVec = @(A) reshape(A, [], 1);
            
            v = zeros(numel(xx)+4, 3);
            
            v(1:numel(xx),:) = [ colVec(xx) colVec(yy) colVec(zz) ];
            
            zBot = obj.zBotFunc(params);
            
            % bottom vertices
            v(end-3:end,:) = [xx(1,1) yy(1,1) zBot;
                xx(end,1) yy(end,1) zBot;
                xx(1,end) yy(1,end) zBot;
                xx(end,end) yy(end,end) zBot];
            
            % reshape into col vec
            v = colVec(v');
        end
        
        function f = faces(obj, v, nx, ny)
            
            x = v(1:3:end);
            y = v(2:3:end);
            z = v(3:3:end);
            
            numTopVerts = nx*ny;
            
            
            %  faces
            f = zeros(2*(nx-1)*(ny-1), 3);
            
            fCount = 0;
            yStart = 1:nx:numTopVerts;
            
            for yi = 1:ny-1
                for xi = 1:nx-1
                    fCount = fCount + 1;
                    f(fCount,:) = [yStart(yi)+xi-1 yStart(yi)+xi yStart(yi+1)+xi];
                    fCount = fCount + 1;
                    f(fCount,:) = [yStart(yi)+xi-1 yStart(yi+1)+xi yStart(yi+1)+xi-1];
                end
            end

            % bottom faces
            f(fCount+1:fCount+2,:) = numTopVerts + [1 4 2; 1 3 4];
            fCount = fCount + 3;
            
            % Side faces
            % Get indices of four lower verts.
            iBot = numTopVerts + [1 2 3 4];
            
            % Indices of vertices along each side of the top surface.
            % Always along ascending x or y direction, so I need to be
            % careful how I flip them around.
            iLowY = 1:nx;
            iHighY = iLowY + nx*(ny-1);
            
            iLowX = 1:nx:nx*ny;
            iHighX = iLowX + nx - 1;
            
            % Side faces!  Sort appropriately to get outward-facing tris.
            
            iSideLowX = [iLowX iBot([3 1])];
            iSideHighX = [iHighX(end:-1:1) iBot([2 4])];
            iSideLowY = [iLowY(end:-1:1) iBot([1 2])];
            iSideHighY = [iHighY iBot([4 3])];
            
            fLowX = iSideLowX(obj.triangulate([-y(iSideLowX), z(iSideLowX)]));
            fHighX = iSideHighX(obj.triangulate([y(iSideHighX), z(iSideHighX)]));
            fLowY = iSideLowY(obj.triangulate([x(iSideLowY), z(iSideLowY)]));
            fHighY = iSideHighY(obj.triangulate([-x(iSideHighY), z(iSideHighY)]));
            
            f = [f; fLowX; fHighX; fLowY; fHighY];
        end
        
        function f = triangulate(obj, verts)
            assert(size(verts, 2) == 2);
            profile = verts;
            
            numVerts = size(profile, 1);
            outerConstraint = [(1:numVerts)', [2:numVerts, 1]'];
            dt = DelaunayTri(profile, outerConstraint);
            inside = inOutStatus(dt);
            
            f = dt.Triangulation(inside,:);
        end
        
        function m = meshes(obj, varargin)
            
            import t6.model.*
            
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            [v nx ny] = obj.vertices(params);
            f = obj.faces(v, nx, ny);
            
            vFunc = @(p) obj.vertices(p);
            myJacobian = jacobian(vFunc, params);
            
            m = { Mesh(v, f, myJacobian, obj.permittivity, ...
                obj.permeability) };
        end
        
        
    end
    
end
