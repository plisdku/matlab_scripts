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
                
                if strcmpi(instruction, 'Simplify')
                    % Use: 'Simplify', costFn
                    % Use: 'Simplify'

                    if iCur+1 <= numInstructions && ...
                       ~ischar(obj.instructions{iCur+1})
                        costFn = obj.instructions{iCur+1};
                        iCur = iCur + 2;

                        [VV, vertices, T, key, crease{:}] = ...
                            VVMesh.simplify(VV, vertices, costFn, crease{:});
                    else
                        [VV, vertices, T, key, crease{:}] = ...
                            VVMesh.simplify(VV, vertices, [], crease{:});
                        iCur = iCur + 1;
                    end
                    
                    if ~isempty(refineThese)
                        refineThese = key(refineThese);
                    end
                    
                elseif strcmpi(instruction, 'Subdivide')
                    % Use: 'Subdivide'
                    % Use: 'Subdivide', vertexPredicate
                    
                    %fprintf('Subdivide (#%i)\n', iCur);
                    if iCur+1 <= numInstructions && ...
                        ~ischar(obj.instructions{iCur+1})
                        predicate = obj.instructions{iCur+1};
                        iCur = iCur + 2;
                        
                        flags = false(size(vertices,1),1);
                        for vv = 1:size(vertices,1)
                            flags(vv) = predicate(vertices(vv,:));
                        end
                        refineThese = find(flags);
                    else
                        iCur = iCur + 1;
                    end
                    
                    if isempty(refineThese)
                        [VV, vertices, T, crease{:}] = ...
                            VVMesh.loopSubdivision(VV, vertices, crease{:});
                    else
                        [VV, vertices, T, refineThese, crease{:}] = ...
                            VVMesh.loopSubdivisionAdaptive(VV, vertices,...
                            refineThese, crease{:});
                    end
                    
                end
                
                if isempty(T_total)
                    T_total = T;
                else
                    T_total = T*T_total;
                end
                
            end
            
            % Transform the Jacobian!
            j0 = mChild{1}.jacobian;
            
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