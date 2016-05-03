function [v2, f2] = reduceVertexList(v1, f1)
% reduceVertexList    Cull extraneous vertices from a face-vertex mesh
%
% [v2,f2] = reduceVertxeList(v1,f1) returns vertices v2 from v1 that are
% referenced in f1, and re-indexes the faces f1 as f2, referring to rows of
% v2.

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

