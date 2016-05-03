function [allVertices, allFaces, totalJacobian] = groupMeshes(varargin)

numMeshes = nargin;
meshes = varargin;

%%
%meshes = {m1, m2};
%numMeshes = length(meshes);
%%
totalNumVertices = sum(cellfun(@(m) m.numVertices, meshes));
totalNumFaces = sum(cellfun(@(m) m.numFaces, meshes));

firstVertexIndex = cumsum([1, cellfun(@(m) m.numVertices, meshes)]);

%% Group the vertices

allVertices = zeros(totalNumVertices, 3);

iNext = 1;
for ii = 1:length(meshes)
    allVertices(iNext:(iNext+meshes{ii}.numVertices-1),:) = meshes{ii}.patchVertices;
    iNext = iNext + meshes{ii}.numVertices;
end

%% Group the faces

allFaces = zeros(totalNumFaces, 3);

iNext = 1;
for ii = 1:length(meshes)
    allFaces(iNext:(iNext+meshes{ii}.numFaces-1),:) = meshes{ii}.faces ...
        + (firstVertexIndex(ii)-1);
    iNext = iNext + meshes{ii}.numFaces;
end

%% Group the jacobians

numParams = max(cellfun(@(m) size(m.jacobian, 2), meshes));

totalJacobian = sparse(3*totalNumVertices, numParams);

for ii = 1:length(meshes)
    [row,col,val] = find(meshes{ii}.jacobian);
    totalJacobian = totalJacobian + sparse(row + 3*(firstVertexIndex(ii)-1),...
        col, val, size(totalJacobian,1), size(totalJacobian,2));
end

