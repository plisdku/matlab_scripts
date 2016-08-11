function [outVertices, outFaces, outParents] = trianglePaternity3d(...
    vChild, fChild, vParents, fParents)
% trianglePaternity3d   Identify intersections of a triangle with several parent triangles
%
% [v, f, parents] = trianglePaternity3d(vChild, fChild, vParents, fParents) returns a
% triangulation of the child triangle, labeled by the index of intersecting parent
% triangles.  The parent and child triangles are assumed to all be coplanar, and the
% analysis is performed in a projection plane.
%
% Parts of the child triangle that intersect no parent triangle are excluded from the
% returned triangulation, so it's up to the caller to make sure the child is really
% covered by the parents.

%% Make a basis matrix.

v1 = vChild(2,:) - vChild(1,:);
v2 = vChild(3,:) - vChild(1,:);
nv = cross(v1, v2);

basis = [v1', v2', nv'];

%% Project down to 2d
% We go into the basis of the child triangle then discard the z-component,
% which should be the same for all the triangles anyway.

% Put origin at vChild(1,:)
vChild_tx = bsxfun(@plus, vChild, -vChild(1,:));
vParent_tx = bsxfun(@plus, vParents, -vChild(1,:));

vChild_transformed = transpose(basis \ vChild_tx');
vParents_transformed = transpose(basis \ vParent_tx');

%%

[v2d, outFaces, outParents] = ll.inherit.trianglePaternity(...
    vChild_transformed(:,1:2), fChild, vParents_transformed(:,1:2), fParents);

%% Inflate back to 3d in the original coordinate system

v3d = [v2d, zeros(size(v2d,1),1)];
outVertices = bsxfun(@plus, transpose(basis * v3d'), vChild(1,:));


