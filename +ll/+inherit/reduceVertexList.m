function [v2, f2] = reduceVertexList(v1, f1)
% reduceVertexList    Cull extraneous vertices from a face-vertex mesh
%
% [v2,f2] = reduceVertexList(v1,f1) creates a mesh equivalent to the input mesh,
% guaranteed to have no unused vertices.

%%
%v1 = allVertices;
%f1 = allFaces(10:12,:);
%%

%[usedVertexIndices, iIntoOld, iIntoNew] = unique(f1(:));
[usedVertexIndices, ~, iIntoNew] = unique(f1(:));
% now usedVertexIndices = f1(:)(iIntoOld)
% and f1(:) = usedVertexIndices(iIntoNew)

f2 = reshape(iIntoNew, size(f1));
v2 = v1(usedVertexIndices,:);

