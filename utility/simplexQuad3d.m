function [xx ww nv bc] = simplexQuad3d(N, vertexRows)
% [xx ww normalVector baryCoords] = simplexQuad3D(Npoints, verticesAsRows)

%row = @(A) reshape(A, 1, []);
%xIn = row(xIn);
%yIn = row(yIn);
%zIn = row(zIn);
%vertsIn = [xIn; yIn; zIn];

vertsIn = vertexRows';

%% Calculate the normal vector to the triangle

v1 = vertsIn(:,2) - vertsIn(:,1);
v2 = vertsIn(:,3) - vertsIn(:,1);
nv = cross(v1, v2);
nv = nv/norm(nv);

columnBasis = [null(nv') nv];

% this is a projection of the vertices, with the first vertex at (0,0).
xyz = columnBasis \ (vertsIn - repmat(vertsIn(:,1), [1,3]));

%% Get quadrature points

[xx2d ww] = simplexquad(N, xyz(1:2,:)');
ww = ww';

%% Convert to barycentric coordinates
% That is, decompose xx2d into a sum of vert(2) and vert(3) in the plane,
% and then vert(1) takes the rest.

bcYZ = xyz(1:2, 2:3) \ xx2d';
bc = [1-sum(bcYZ,1); bcYZ];

%% The x and y in the plane...

nodes2d = xx2d(:,1:2)';

%x = xx2d(:,1)';
%y = xx2d(:,2)';

xx = bsxfun(@plus, (columnBasis(:,1:2) * nodes2d), vertsIn(:,1));
