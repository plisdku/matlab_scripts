function [allSubVertices, allSubFaces, allSubJacX, allSubJacY, allSubJacZ] ...
    = inheritJacobian(child, varargin)
% inheritJacobian   Propagate mesh Jacobian matrices through the result of
% Boolean operations
%
% [v, f, Dvx, Dvy, Dvz] = inheritJacobian(child, parent1, parent2, ...)
%
% where child, parent1, parent2 etc. are t6.model.Mesh objects.
%
% Then v (Nx3) and f (Mx3) give vertices and faces on the surface of the
% child mesh, with the Jacobian of each vertex in Dvx, Dvy and Dvz.
% The faces of the child mesh may be partitioned and/or re-triangulated, 
% and the the list of vertices can be expected to include repeats.

[allVertices, allFaces, allJacobian] = ll.inherit.groupMeshes(varargin{:});
allJacobianX = allJacobian(1:3:end,:);
allJacobianY = allJacobian(2:3:end,:);
allJacobianZ = allJacobian(3:3:end,:);

% Determine which faces of the child are coplanar and intersecting with
% faces of the parent meshes.
[iNew, iOld] = neflab.paternity(allVertices, allFaces, child.patchVertices, child.faces);
oldNew = sparse(iOld, iNew, ones(size(iNew)));

%%

% Cache the Nx3 child vertex list, which is recalculated on each call
% to patchVertices().
childVertices = child.patchVertices();

% Let's assume that for N child faces I will get no more than 3*N
% sub-faces.

numExpectedSubFaces = 4*size(child.faces,1);
allSubVertices = zeros(3*numExpectedSubFaces, 3);
allSubFaces = zeros(numExpectedSubFaces,3);
allSubParents = zeros(numExpectedSubFaces,1);
iNextSubVert = 1;
iNextSubFace = 1;

numChildFaces = size(child.faces,1);
for iNew = 1:numChildFaces
    
    iOldFaces = find(oldNew(:,iNew));
    
%     figure(2); clf
%     patch('Faces', child.faces, 'Vertices', child.patchVertices, ...
%         'FaceColor', 'b', 'FaceAlpha', 1)
%     view(3); axis image; camlight right; title('Union');
%     xlabel('x'); ylabel('y'); zlabel('z')
%     hold on
% 
%     ax = axis;
% 
%     figure(3); clf
%     patch('Faces', allFaces(iOldFaces,:), 'Vertices', allVertices, ... 
%        'FaceColor', 'r', 'FaceAlpha', 0.25);
%     patch('Faces', child.faces(iNew,:), 'Vertices', child.patchVertices, ...
%        'FaceColor', 'b', 'FaceAlpha', 0.25);
%     view(3); axis image; camlight right; title('Union');
%     xlabel('x'); ylabel('y'); zlabel('z')
%     axis(ax);
    
    % Isolate the relevant parent faces.  This should help avoid
    % unnecessary edge splittings, but I don't think the step is strictly
    % necessary.
    [vOldReduced, fOldReduced] = ll.inherit.reduceVertexList(allVertices, allFaces(iOldFaces,:));
    [vNewReduced, fNewReduced] = ll.inherit.reduceVertexList(childVertices, child.faces(iNew,:));
    
    [subVertices, subFaces, subParentsReduced] = ll.inherit.trianglePaternity3d(...
        vNewReduced, fNewReduced, vOldReduced, fOldReduced);
    subParents = iOldFaces(subParentsReduced);
    
    
%     % Plot 'em!
%     figure(2);
%     hold on
%     colors = autumn(size(subFaces,1)+2);
%     colors = colors(2:end-1,:);
%     for tt = 1:size(subFaces,1)
%         v = subVertices(subFaces(tt,:),:);
%         
%         % Inset v a bit
%         centroid = mean(v);
%         v = bsxfun(@plus, 0.9*v, 0.1*centroid);
%         
%         plot3(v([1 2 3 1],1), v([1 2 3 1], 2), v([1 2 3 1], 3), ...
%             'Color', colors(tt,:), 'linewidth', 4);
%     end
    
    iNextNextSubVert = iNextSubVert + size(subVertices,1);
    iNextNextSubFace = iNextSubFace + size(subFaces,1);
    
    allSubVertices(iNextSubVert:(iNextNextSubVert-1),:) = ...
        subVertices;
    allSubFaces(iNextSubFace:(iNextNextSubFace-1),:) = subFaces + (iNextSubVert-1);
    allSubParents(iNextSubFace:(iNextNextSubFace-1)) = subParents;
    
    iNextSubVert = iNextNextSubVert;
    iNextSubFace = iNextNextSubFace;
    
    if iNextSubVert > size(allSubVertices,1)
        warning('Did not allocate enough sub vertex space');
    end
    
    if iNextSubFace > size(allSubFaces,1)
        warning('Did not allocate enough sub face space');
    end
    
    %pause
end

allSubVertices = allSubVertices(1:(iNextSubVert-1),:);
allSubFaces = allSubFaces(1:(iNextSubFace-1),:);
allSubParents = allSubParents(1:(iNextSubFace-1),:);

%% Cull the redundant vertices.
% If I improve the code above that partitions the child faces, I can
% perhaps avoid redundant vertices.

[allSubVertices, ~, iIntoUniqueVerts] = unique(allSubVertices, 'rows');
iUsedVertices = allSubFaces(:);
allSubFaces = reshape(iIntoUniqueVerts(iUsedVertices), size(allSubFaces));

%%

numSubFaces = size(allSubFaces,1);
numSubVertices = size(allSubVertices,1);
numParameters = size(allJacobian,2);

% Make some room to build the sparse Jacobians.
% Pre-allocating should hopefully speed up the gross matrix-building
% I'm doing on a per-vertex basis.
estimatedNonzeros = round(1.5 * nnz(allJacobian) * ...
    (numSubVertices / size(allVertices,1)));
allSubJacX = spalloc(numSubVertices, numParameters, estimatedNonzeros);
allSubJacY = spalloc(numSubVertices, numParameters, estimatedNonzeros);
allSubJacZ = spalloc(numSubVertices, numParameters, estimatedNonzeros);

for ss = 1:numSubFaces
    
    % Put the subtriangle in child-triangle barycentric coordinates.
    % Find the Jacobian for each new vertex.
    % 
    
    iSubVerts = allSubFaces(ss,:);
    v = allSubVertices(iSubVerts,:);
    
    parentVertexIndices = allFaces(allSubParents(ss),:);
    parentVertices = allVertices(parentVertexIndices,:);
    parentJacX = allJacobianX(parentVertexIndices,:);
    parentJacY = allJacobianY(parentVertexIndices,:);
    parentJacZ = allJacobianZ(parentVertexIndices,:);
    
    [subJacX, subJacY, subJacZ] = ll.inherit.subTriangleJacobian(parentVertices,...
        parentJacX, parentJacY, parentJacZ, ...
        v);
    
    allSubJacX(iSubVerts, :) = subJacX;
    allSubJacY(iSubVerts, :) = subJacY;
    allSubJacZ(iSubVerts, :) = subJacZ;
    
end
