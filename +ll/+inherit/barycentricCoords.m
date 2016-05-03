function baryCoords = barycentricCoords( basisVertices, points )
% barycentricCoords   Express points on a triangle in 3D space as a
% weighted sum of the triangle's corner points.
%
% baryCoords = barycentricCoords(basisVertices, points)
%
% takes the triangle vertices as three columns of basisVertices (3x3), and
% the points to decompose in a 3xN array.

A = [basisVertices(:,1) - basisVertices(:,3), ...
    basisVertices(:,2) - basisVertices(:,3)];

coeffs2d = A \ (points - repmat(basisVertices(:,3), 1, size(points,2)));

baryCoords = [coeffs2d; 1 - coeffs2d(1,:) - coeffs2d(2,:)];

