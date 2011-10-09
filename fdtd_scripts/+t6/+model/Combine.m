classdef Combine < t6.model.Node
    
    properties
        children = {};
    end
    
    
    methods
        function obj = Combine(varargin)
            if nargin > 0
                obj.children = varargin;
            end
        end
        
        function m = meshes(obj, varargin)
            import t6.model.*;
            
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            allBounds = [];
            
            % 1.  Determine total number of vertices and faces
            doOverlap = @(r1,r2) all(r1(4:6) > r2(1:3) & r1(1:3) < r2(4:6));
            
            numVerts = 0;
            numFaces = 0;
            
            ourPermittivity = '';
            ourPermeability = '';
            
            for cc = 1:length(obj.children)
                childMeshes = obj.children{cc}.meshes(params);
                
                for mm = 1:length(childMeshes)
                    allBounds(end+1,:) = childMeshes{mm}.bounds();
                    numVerts = numVerts + length(childMeshes{mm}.vertices)/3;
                    numFaces = numFaces + length(childMeshes{mm}.faces);
                    
                    ourPermittivity = childMeshes{mm}.permittivity;
                    ourPermeability = childMeshes{mm}.permeability;
                end
            end
            
            % 2.  Put up a warning if any child meshes have overlapping
            % bounding boxes
            
            for mm = 1:size(allBounds,1)
                for nn = mm+1:size(allBounds,1)
                    if doOverlap(allBounds(mm,:), allBounds(nn,:))
                        warning(['Some combined meshes (%i and %i) ', ...
                            'have overlapping bounding boxes'], mm, nn);
                    end
                end
            end
            
            % 3.  Report an error if the child meshes do not all have the
            % same permittivity and permeability.
            
            for cc = 1:length(obj.children)
                childMeshes = obj.children{cc}.meshes(params);
                
                for mm = 1:length(childMeshes)
                    if ~strcmp(childMeshes{mm}.permittivity, ourPermittivity)
                        error('Not all child meshes have same permittivity');
                    end
                    
                    if ~strcmp(childMeshes{mm}.permeability, ourPermeability)
                        error('Not all child meshes have same permeability');
                    end
                end
            end
            
            % 4. Concatenate all the vertices, faces and jacobians.
            allVertices = zeros(numVerts*3, 1);
            allFaces = zeros(numFaces,3);
            allJacobians = sparse(numVerts*3, size(params,1));
            
            nextVert = 1;
            nextFace = 1;
            
            for cc = 1:length(obj.children)
                childMeshes = obj.children{cc}.meshes(params);
                
                for mm = 1:length(childMeshes)
                    nv = size(childMeshes{mm}.vertices,1);
                    nf = size(childMeshes{mm}.faces, 1);
                    allVertices(nextVert:nextVert + nv - 1) = ...
                        childMeshes{mm}.vertices;
                    allFaces(nextFace:nextFace + nf - 1, :) = ...
                        childMeshes{mm}.faces + (nextVert - 1)/3;
                    allJacobians(nextVert:nextVert + nv - 1, :) = ...
                        childMeshes{mm}.jacobian;
                    
                    nextVert = nextVert + nv;
                    nextFace = nextFace + nf;
                end
            end
            
            m = { Mesh(allVertices, allFaces, allJacobians, ...
                ourPermittivity, ourPermeability) };
            
        end
    end
end
        