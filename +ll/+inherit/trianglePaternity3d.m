function [outVertices, outFaces, outParents] = trianglePaternity3d(...
    vChild, fChild, vParents, fParents)

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


