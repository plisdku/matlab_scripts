function nv = facetNormal(vertices, face)
% calculate normal vector to a facet...

numVerts = numel(face);
orderedVerts = vertices(face,:);

% Calculate an oriented area vector.

nv = [0 0 0];

for vv = 2:(numVerts-1)
    v1 = orderedVerts(vv,:) - orderedVerts(1,:);
    v2 = orderedVerts(vv+1,:) - orderedVerts(1,:);
    
    nv = nv + cross(v1, v2);
end

nv = nv / norm(nv);
