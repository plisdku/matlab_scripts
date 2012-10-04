classdef ExtrudePath < t6.model.Node
    
    properties
        permittivity = 'none';
        permeability = 'none';
        xFunc = [];
        yFunc = [];
        pathFunc = [];  % map [0 1] to [x y z]
        fwdFunc = [];
        uFunc = [];    % map [0 1] to [ux uy uz]
        vFunc = [];
        t = [];
    end
    
    methods
        function obj = ExtrudePath(varargin)
            
            X.X = [];
            X.Y = [];
            X.Permittivity = '';
            X.Permeability = '';
            X.Path = [];
            X.U = [];
            X.V = [];
            X.T = linspace(0, 1, 10);
            X = parseargs(X, varargin{:});
            
            if isempty(X.U) && isempty(X.V)
                error('Neither U nor V was provided');
            end
            
            if ~isa(X.X, 'function_handle')
                error('X must be a function handle');
            end
            
            if ~isa(X.Y, 'function_handle')
                error('Y must be a function handle');
            end
            
            if ~isa(X.Path, 'function_handle')
                error('Path must be a function handle');
            end
            
            if ~isempty(X.U) && ~isa(X.U, 'function_handle')
                error('U must be a function handle');
            end
            
            if ~isempty(X.V) && ~isa(X.V, 'function_handle')
                error('V must be a function handle');
            end
            
            if ~isa(X.T, 'numeric')
                error('T must be an array of numbers in [0, 1]');
            end
            
            obj.permittivity = X.Permittivity;
            obj.permeability = X.Permeability;
            obj.xFunc = X.X;
            obj.yFunc = X.Y;
            obj.pathFunc = X.Path;
            obj.uFunc = X.U;
            obj.vFunc = X.V;
            obj.t = reshape(X.T, 1, []);
        end
        
        function [xBot yBot tris] = bottomFace(obj, params)
            
            % Make the bottom face
            xBot = obj.xFunc(params);
            yBot = obj.yFunc(params);
            
            if size(xBot, 2) ~= 1 || size(yBot, 2) ~= 1
                error('X and Y must return column vectors');
            end
            
            if nargout > 2

                profile = [xBot yBot];
                outerConstraint = [(1:numel(xBot))', [2:numel(xBot), 1]'];
                dt = DelaunayTri(profile, outerConstraint);
                inside = inOutStatus(dt);

                tris = dt.Triangulation(inside,:);
            end
            
            %DxBot = t6.model.jacobian(obj.xFunc, params);
            %DyBot = t6.model.jacobian(obj.yFunc, params);
            
        end
        
        function allFaces = faces(obj, tris, numEdges)
            
            v0 = 1:numEdges;
            v1 = [2:numEdges, 1];
            v2 = v1 + numEdges;
            v3 = v0 + numEdges;

            ringFaces = [v0' v1' v2'; v0' v2' v3'];

            numEndFaces = size(tris, 1);
            numRings = length(obj.t)-1;
            allFaces = zeros(numEdges*2*numRings + 2*numEndFaces, 3);

            facesPerRing = 2*numEdges;
            for rr = 1:numRings
                allFaces((rr-1)*facesPerRing + (1:facesPerRing), :) = ...
                    ringFaces + (rr-1)*numEdges;
            end

            % Add end faces.

            allFaces(numRings*facesPerRing + (1:numEndFaces), :) = fliplr(tris);
            allFaces(numRings*facesPerRing + numEndFaces + (1:numEndFaces), :) = ...
                tris + numRings*numEdges;
            
        end
        
        function v = allVertices(obj, params)
            
            [xBot yBot] = obj.bottomFace(params);
            v = obj.vertices(params, xBot, yBot);
            
        end
        
        function v = vertices(obj, params, xBot, yBot)
            
            % 1. calculate profile vertices
            % 2. calculate spine path
            % 3. calculate the rotation matrix at each position
            % 4. arrange vertices as single column (right?)
            
            path = cell2mat(arrayfun(@(t) obj.pathFunc(params,t), obj.t, ...
                'UniformOutput', false));
            fwd = centeredDiff(path, 2);
            fwd = bsxfun(@times, fwd, 1./sqrt(sum(fwd.^2,1)));
            
            if ~isempty(obj.uFunc)
                u = cell2mat(arrayfun(@(t) obj.uFunc(params,t), obj.t, ...
                    'UniformOutput', false));
            end
            
            if ~isempty(obj.vFunc)
                v = cell2mat(arrayfun(@(t) obj.vFunc(params,t), obj.t, ...
                    'UniformOutput', false));
            end
            
            if ~exist('u', 'var')
                u = cross(v, fwd, 1);
            end
            
            if ~exist('v', 'var')
                v = cross(fwd, u, 1);
            end
            
            % Build local basis func!
            
            uTensor = reshape(u, 3, 1, []);
            vTensor = reshape(v, 3, 1, []);
            pTensor = reshape(path, 3, 1, []);
            
            A = cat(2, uTensor, vTensor, pTensor);
            xy = [xBot yBot ones(size(xBot))];
            
            stretchOut = @(A) A(:);
            v = stretchOut(t6.multTensor(A, xy, 2));
        end
        
        
        function m = meshes(obj, varargin)
            import t6.model.*
            
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            [xBot yBot tris] = obj.bottomFace(params);
            
            vertFunc = @(p) obj.allVertices(p);
            
            myVerts = vertFunc(params);
            myFaces = obj.faces(tris, length(xBot));
            
            myJacobian = jacobian(vertFunc, params);
            
            m = { Mesh(myVerts, myFaces, myJacobian, obj.permittivity, ...
                obj.permeability) };
        end
        
        
        
    end % methods
    
    
end
            
            
            
            