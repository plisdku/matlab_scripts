classdef ProcessMesh < t6.model.Node
    
    properties
        childNode = [];
        instructions = {};
        creases = {};
        refineVertices = [];
        permittivity = 'none';
        permeability = 'none';
    end
    
    
    methods
        
        % ProcessMesh(childMesh, creases, refineVertices)
        % ProcessMesh(childMesh, creases, refineVertices, ...
        %   'Simplify', costFn, 'Subdivide', 'Subdivide', ...
        %   'Simplify', ...)
        function obj = ProcessMesh(child, creases, refineVertices, varargin)
            obj.childNode = child;
            obj.creases = creases;
            obj.refineVertices = refineVertices;
            obj.instructions = varargin;
            obj.permittivity = child.permittivity;
            obj.permeability = child.permeability;
        end
        
        function [m, T_total] = meshes(obj, varargin)
            
            if nargin > 1
                params = varargin{1};
            else
                params = [];
            end
            
            mChild = obj.childNode.meshes(params);
            
            vertices = mChild{1}.patchVertices;
            VV = VVMesh.fv2vv(mChild{1}.faces);
            crease = obj.creases;
            refineThese = obj.refineVertices;
            
            % Apply simplifications:
            iCur = 1;
            
            numInstructions = numel(obj.instructions);
            T_total = [];
            while iCur <= numInstructions
                
                instruction = obj.instructions{iCur};
                
                if iCur+1 <= numInstructions && ...
                    ~ischar(obj.instructions{iCur+1})
                    
                    argument = obj.instructions{iCur+1};
                    numArguments = 1;
                    iCur = iCur + 2;
                else
                    numArguments = 0;
                    iCur = iCur + 1;
                end
                
                fprintf('Handling %s\n', instruction);
                
                if strcmpi(instruction, 'Simplify')
                    % Use: 'Simplify', costFn
                    % Use: 'Simplify'
                    
                    if numArguments == 1;
                        costFn = argument;
                        [VV, vertices, T, key, crease{:}] = ...
                            VVMesh.simplify(VV, vertices, costFn, crease{:});
                    else
                        [VV, vertices, T, key, crease{:}] = ...
                            VVMesh.simplify(VV, vertices, [], crease{:});
                    end
                    
                    if ~isempty(refineThese)
                        refineThese = key(refineThese);
                    end
                    
                elseif strcmpi(instruction, 'Subdivide')
                    % Use: 'Subdivide'
                    % Use: 'Subdivide', vertexPredicate
                    
                    %fprintf('Subdivide (#%i)\n', iCur);
                    if numArguments == 1
                        flags = false(size(vertices,1),1);
                        for vv = 1:size(vertices,1)
                            flags(vv) = argument(vertices(vv,:));
                        end
                        refineThese = find(flags);
                        
                        fprintf('There are %i verts to subdivide.\n', ...
                            numel(refineThese));
                    end
                    
                    if isempty(refineThese)
                        [VV, vertices, T, crease{:}] = ...
                            VVMesh.loopSubdivision(VV, vertices, crease{:});
                    else
                        [VV, vertices, T, refineThese, crease{:}] = ...
                            VVMesh.loopSubdivisionAdaptive(VV, vertices,...
                            refineThese, crease{:});
                    end
                    
                elseif strcmpi(instruction, 'Refine')
                    % Use: 'Refine'
                    % Use: 'Refine', vertexPredicate
                    
                    if numArguments == 1
                        flags = false(size(vertices,1),1);
                        for vv = 1:size(vertices,1)
                            flags(vv) = argument(vertices(vv,:));
                        end
                        refineThese = find(flags);
                        
                        fprintf('There are %i verts to refine.\n', ...
                            numel(refineThese));
                    end
                    
                    if numArguments == 0
                        [VV2, vertices, ~, perturbFlags, T] = ...
                            VVMesh.loopRefine(VV, vertices);
                    elseif ~isempty(refineThese)
                        [VV2, vertices, ~, perturbFlags, T] = ...
                            VVMesh.loopRefine(VV, vertices, refineThese);
                    else
                        VV2 = VV;
                        T = speye(size(VV,1));
                    end
                    
                    % Update crease vertex indices so the vertices of all
                    % edges that were split are also included.
                    if ~isequal(VV2, VV)
                    for cc = 1:numel(crease)
                        crease{cc} = VVMesh.subdividedCrease(...
                            crease{cc}, perturbFlags(crease{cc}),...
                            VV, VV2);
                        % T is correct regardless.
                    end
                    end
                    
                    VV = VV2;
                    
                    % Expand the refinement region for the next refinement
                    % level.
                    if ~isempty(refineThese)
                        refineThese = union(refineThese, ...
                            VVMesh.neighborhood(refineThese, VV, 1, 2));
                    end
                    
                end % if simplify, subdivide, refine
                
                if isempty(T_total)
                    T_total = T;
                else
                    T_total = T*T_total;
                end
                
            end
            
            
            % Transform the Jacobian!
            j0 = mChild{1}.jacobian;
            
            if isempty(T_total)
                T_total = speye(size(j0,1)/3);
            end
            
            jx = T_total*j0(1:3:end,:);
            jy = T_total*j0(2:3:end,:);
            jz = T_total*j0(3:3:end,:);
            jxyz = vertcat(jx,jy,jz);
            
            % All this horrible indexing makes me think it may have been an
            % error to sort [x0 y0 z0 x1 y1 z1] instead of putting all the
            % xs first, then ys.  Barf.
            rowOrdering = reshape(1:numel(vertices),[],3)';
            myJacobian = jxyz(rowOrdering,:);
            
            % Make the faces
            myFaces = VVMesh.vv2fv(VV);
            
            % Arrange the vertices
            myVertices = zeros(numel(vertices),1);
            myVertices(1:3:end) = vertices(:,1);
            myVertices(2:3:end) = vertices(:,2);
            myVertices(3:3:end) = vertices(:,3);
            
            m = { t6.model.Mesh(myVertices, myFaces, myJacobian, ...
                obj.permittivity, obj.permeability) };
            
        end
        
            
    end
    
end