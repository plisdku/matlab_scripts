function shells = divideIntoShells(faces)
% Given polyhedral faces, divide into a bunch of disjoint shells.
% a simple polyhedron has one shell.
%
% Two faces are adjacent if they share two vertices.
%
% I'll create a vertex adjacency matrix first, such that M(face, vert) = 1
% iff the face has that vertex index in it.

% Make the face-to-vertices adjacency matrix:

numFaces = size(faces,1);

rows_faces = ceil( (1:3*numFaces)/3 );
cols_verts = faces';

f2v = sparse(rows_faces, cols_verts(:), ones(size(rows_faces)));

%f2f = f2v * f2v' - 3*eye(numFaces);
%faceAdjacency = (f2f == 2);

faceAdjacency = (f2v * f2v') >= 2; % this INCLUDES the diags, on purpose

% Now divide into freaking shells!
% I can use the extremely awesome dmperm function to do this!!
[permutedFaceNumbers, ~, whichNodes] = dmperm(faceAdjacency);

numShells = numel(whichNodes)-1;

shells = cell(numShells, 1);

for ss = 1:numShells
    shells{ss} = faces(permutedFaceNumbers(whichNodes(ss):whichNodes(ss+1)-1),:);
end

%% A driver function that I can run cellwise to test things out
function runThisFunction
%%

[shells, M] = divideIntoShells([1 2 3; 3 2 5])

figure(1); clf
spy(M);

figure(2); clf
imagesc(M, [0 2]); colorbar

% Dulmage-Mendelsohn decomposition figures out the connected components.

[p,q,r,s] = dmperm(M);
% p  permuted list of nodes
% r  which nodes of p belong to the same connected component
%    e.g. p(r(1):r(2)-1)







