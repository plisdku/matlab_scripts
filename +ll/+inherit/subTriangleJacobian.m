function [Dwx, Dwy, Dwz] = subTriangleJacobian(v, Dvx, Dvy, Dvz, w)
% subTriangleJacobian   Compute Jacobians for vertices bound to a movable
% triangle.
%
% [Dwx, Dwy, Dwz] = subTriangleJacobian(v, Dvx, Dvy, Dvz, w) will express
% w in barycentric coordinates (weighted combination of points v) and
% return weighted combinations of the Jacobians Dvx, Dvy, and Dvz of the
% parent triangle.
%
% w is 3x3, indexed w(vertexIndex, xyz).
% Dvx is 3xN, indexed Dvx(vertexIndex, parameterIndex).
% Dvy, Dvz are like Dvx.
%
% Dwx, Dwy, Dwz are also like Dvx.

% we do sparse and non-sparse versions just so the returned Jacobians
% are sparse when the input Jacobians are sparse.
if issparse(Dvx)
    baryWeights = sparse(ll.inherit.barycentricCoords(v', w'));
else
    baryWeights = ll.inherit.barycentricCoords(v', w');
end

Dwx = baryWeights' * Dvx;
Dwy = baryWeights' * Dvy;
Dwz = baryWeights' * Dvz;